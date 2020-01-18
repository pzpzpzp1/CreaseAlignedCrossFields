#include "mex.hpp"
#include "mexAdapter.hpp"
#include <cstdlib>
#include "fusion.h"
#include <csignal>
#include <ctime>

/*
REMEMBER MOSEK'S RESHAPE IS ROW-MAJOR WHICH IS DIFFERENT FROM MATLAB'S COLUMN-MAJOR.
*/

using namespace monty;
using namespace mosek::fusion;
using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
    ArrayFactory factory;
	Array weights;
	std::vector<double> timings;
	char buffer[300];
	Model::t M;
	std::clock_t start;
public:

	void tic() {
		start = std::clock();
	}

	double toc(void) {
		return toc(start);
	}

	double toc(clock_t base) {
		auto duration = (std::clock() - base) / (double)CLOCKS_PER_SEC;

		sprintf(buffer, "Elapsed time: %f seconds\n", duration);
		displayCharArrayOnMATLAB(buffer);

		return duration;
	}

	TypedArray<double> solveSoftCrossFields(Matrix::t& A, Matrix::t& D, Matrix::t& b, Matrix::t& a, int p, int nEdges, int nFaces, double n, Matrix::t& c, int nTotalFaces)
	{
		// Create Variables
		displayStringOnMATLAB("Creating Variables: ");
		M = new Model();
		M->setSolverParam("intpntSolveForm", "dual");
		M->setSolverParam("numThreads", 0);
		auto x = M->variable(9 * nFaces, Domain::unbounded());
		auto y = M->variable(9 * nEdges, Domain::unbounded());
		auto ynorms = M->variable(nEdges, Domain::greaterThan(0.));
		auto z = M->variable(1, Domain::greaterThan(0.));
		timings.push_back(toc()); tic();

		// Objective
		displayStringOnMATLAB("Building Objective: ");
		M->objective(ObjectiveSense::Minimize, z);
		timings.push_back(toc()); tic();

		// Normal alignment constraint
		displayStringOnMATLAB("Normal Alignment. ");		
		if (n == 0) {
			// Ax == b
			auto Axmb = Expr::sub(Expr::mul(A, x), b);
			displayStringOnMATLAB("N=0. ");
			Constraint::t normalAlignment = M->constraint(Axmb, Domain::equalsTo(0));
		}
		else {
			// reshape(Ax,7,[]) - repmat(b,1,nFaces) \in \Qcone
			
			auto Ax = Expr::mul(A, x);
			auto Axmb = Expr::sub(Ax, b);
			auto AxmbReshaped = Expr::reshape(Axmb, nFaces, 7);

			displayStringOnMATLAB("N~=0. ");
			auto nRepeated = Expr::repeat(Expr::constTerm(n),nFaces,0);
			auto normalAlignmentConePoint = Expr::hstack(nRepeated, AxmbReshaped);
			Constraint::t normalAlignment = M->constraint(normalAlignmentConePoint, Domain::inQCone(nFaces, 8));
		}

		// Constraint defining y = Dx+c. 
		displayStringOnMATLAB("Defining y=Dx+c. ");
		Constraint::t defineY = M->constraint(Expr::sub(y, Expr::add(Expr::mul(D, x), c)), Domain::equalsTo(0)); // y=Dx+c

		// Constraint defining ynorms = norms(reshape(y))
		displayStringOnMATLAB("Defining ynorms. ");
		auto multiYnormConstraints = Expr::hstack(ynorms, y->reshape(nEdges, 9));
		auto y2normConstraint = M->constraint(multiYnormConstraints, Domain::inQCone(nEdges, 10));
		timings.push_back(toc()); tic();
		
		// P-norm objective as a cone constraint
		displayStringOnMATLAB("Building pnorm ");
		if (p==1){
			//l1 norm. //z > sum(ynorms.*a)
			displayStringOnMATLAB("P=1: ");
			Constraint::t pnorm = M->constraint(Expr::sub(z, Expr::dot(a,ynorms)), Domain::greaterThan(0)); 
		}
		else if (p==2){
			//l2 norm. use inQCone. z > sqrt(sum(a * ynorms^2))
			displayStringOnMATLAB("P=2: ");
			Matrix::t sqrta = matlabVectorArrayToDenseMatrixSqrt(weights);
			Constraint::t pnorm = M->constraint(Expr::vstack(z, Expr::reshape(Expr::mulElm(ynorms, sqrta), nEdges)), Domain::inQCone()); 
		}
		else if (p==-1){
			// infinity. max. z > max(ynorms). weights have no effect.
			displayStringOnMATLAB("P=inf: ");
			Constraint::t pnorm = M->constraint(Expr::sub(Var::repeat(z, nEdges), ynorms), Domain::greaterThan(0.));
		}
		else {
			// use pnorm cone
			displayStringOnMATLAB("P=Other: ");
			Matrix::t powa = matlabVectorArrayToDenseMatrixPow(weights, p);
			pnorm(M, z, Expr::mulElm(ynorms, powa), p); 
		}
		timings.push_back(toc()); tic();
		
		// Solve
		displayStringOnMATLAB("Solving. ");
		M->solve();
		timings.push_back(toc()); tic();
		
		// Write output
		displayStringOnMATLAB("Writing output. ");
		TypedArray<double> outN = factory.createArray<double>(
			{ static_cast<size_t>(nFaces * 9), static_cast<size_t>(1) });
		try {
			auto result = *(x->level());
			for (int i = 0; i < 9 * nFaces; i++)
			{
				outN[i] = result[i];
			}
			timings.push_back(toc()); tic();
		}
		catch (mosek::fusion::FusionException mosekexcept) {
			disposeM();
			throw(mosekexcept);
		}

		// Dispose Model
		disposeM();

		return outN;
	}

	void disposeM(void) {
		tic();
		displayStringOnMATLAB("Disposing Model.");
		M->dispose();
		timings.push_back(toc()); tic();
	}

	Matrix::t matlabSparseArrayToSparseMatrix(SparseArray<double>& a)
	{
		ArrayDimensions Adim = a.getDimensions();
		int nnz = a.getNumberOfNonZeroElements();
		auto subi = std::shared_ptr<ndarray<int, 1> >(new ndarray<int, 1>(shape(nnz)));
		auto subj = std::shared_ptr<ndarray<int, 1> >(new ndarray<int, 1>(shape(nnz)));
		auto val = std::shared_ptr<ndarray<double, 1> >(new ndarray<double, 1>(shape(nnz)));
		
		int counter = 0;
		for(auto it = a.begin(); it!=a.end(); it++)
		{
			(*val)(counter) = it[0];
			(*subi)(counter) = a.getIndex(it).first;
			(*subj)(counter) = a.getIndex(it).second;
			counter ++;

			if (counter % 30000 == 0) {
				sprintf(buffer, "%d nonzeros processed (%f). Val:%f\n", counter, ((double)counter)/nnz, it[0]);
				displayCharArrayOnMATLAB(buffer);
				if (counter == -1) {
					assert(false);
				}
			}
		}
		
		Matrix::t b = Matrix::sparse(Adim[0], Adim[1], subi, subj, val);
		return b;
	}

	Matrix::t matlabTripletsToSparseMatrix(Array& Di, Array& Dj, Array& Dv, int D0, int D1) {
		ArrayDimensions dims = Di.getDimensions();
		int nnz = dims[0];
		auto subi = std::shared_ptr<ndarray<int, 1> >(new ndarray<int, 1>(shape(nnz)));
		auto subj = std::shared_ptr<ndarray<int, 1> >(new ndarray<int, 1>(shape(nnz)));
		auto val = std::shared_ptr<ndarray<double, 1> >(new ndarray<double, 1>(shape(nnz)));

		for (int i = 0; i < nnz; i++) {
			(*val)(i) = Dv[i];
			(*subi)(i) = Di[i];
			(*subj)(i) = Dj[i];
		}
		Matrix::t b = Matrix::sparse(D0, D1, subi, subj, val);
		return b;
	}

	Matrix::t matlabVectorArrayToDenseMatrix(Array& a)
	{
		ArrayDimensions Adim = a.getDimensions();
		auto val = std::shared_ptr<ndarray<double, 1> >(new ndarray<double, 1>(shape(Adim[0])));
		for(int i = 0; i < Adim[0]; i++)
		{
			(*val)[i] = a[i];
		}
	
		Matrix::t b = Matrix::dense(Adim[0], Adim[1], val);
		return b;
	}
	
	Matrix::t matlabVectorArrayToDenseMatrixSqrt(Array& a)
	{
		ArrayDimensions Adim = a.getDimensions();
		auto val = std::shared_ptr<ndarray<double, 1> >(new ndarray<double, 1>(shape(Adim[0])));
		for(int i = 0; i < Adim[0]; i++)
		{
			(*val)[i] = std::sqrt((double)a[i]);
		}
	
		Matrix::t b = Matrix::dense(Adim[0], Adim[1], val);
		return b;
	}
	
	// a^(1/p)
	Matrix::t matlabVectorArrayToDenseMatrixPow(Array& a, int p)
	{
		ArrayDimensions Adim = a.getDimensions();
		auto val = std::shared_ptr<ndarray<double, 1> >(new ndarray<double, 1>(shape(Adim[0])));
		for(int i = 0; i < Adim[0]; i++)
		{
			(*val)[i] = std::pow((double)a[i], 1.0/p);
		}
	
		Matrix::t b = Matrix::dense(Adim[0], Adim[1], val);
		return b;
	}
	
	// p-norm, p>1
	// t >= \|x\|_p (where p>1), x is a vector variable
	void pnorm(Model::t M, Expression::t t, Expression::t x, double p) {
		int n = (int)x->getSize();
		Variable::t r = M->variable(n);
		M->constraint(Expr::sub(t, Expr::sum(r)), Domain::equalsTo(0.0));
		auto rept = Expr::reshape(Expr::repeat(t, n, 1), n);

		M->constraint(Expr::hstack(rept, r, x), Domain::inPPowerCone(1.0 - 1.0 / p));
	}

	// Start of mex function
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
		// Extract inputs and verify
		timings.clear(); tic();
		Array Ai = std::move(inputs[0]);
		Array Aj = std::move(inputs[1]);
		Array Av = std::move(inputs[2]);
		Array A0 = std::move(inputs[3]);
		Array A1 = std::move(inputs[4]);
		Array b = std::move(inputs[5]);
		Array Di = std::move(inputs[6]);
		Array Dj = std::move(inputs[7]);
		Array Dv = std::move(inputs[8]); 
		Array D0 = std::move(inputs[9]);
		Array D1 = std::move(inputs[10]);
		Array a = std::move(inputs[11]);
		Array P = std::move(inputs[12]);
		Array N = std::move(inputs[13]);
		Array c = std::move(inputs[14]);
		int a0 = A0[0];
		int a1 = A1[0];
		int d0 = D0[0];
		int d1 = D1[0];
		int p = P[0];
		double n = N[0];
		weights = a;
		
        ArrayDimensions bdim = b.getDimensions();
		ArrayDimensions adim = a.getDimensions();
		
		int nEdges = d0/9;
		int nFaces = a1/9;
		int nTotalFaces = a0 / 7;
		assert(7*nFaces==bdim[0]);
		assert(p>=1 || p==-1);
		assert(n>=0);
		assert(nEdges == adim[0]);
		
		displayStringOnMATLAB("Converting Array to Matrix: ");
		Matrix::t D_matrix = matlabTripletsToSparseMatrix(Di, Dj, Dv, d0, d1);
		Matrix::t A_matrix = matlabTripletsToSparseMatrix(Ai, Aj, Av,a0,a1);
		Matrix::t a_matrix = matlabVectorArrayToDenseMatrix(a);
		Matrix::t b_matrix = matlabVectorArrayToDenseMatrix(b);
		Matrix::t c_matrix = matlabVectorArrayToDenseMatrix(c);
		timings.push_back(toc()); tic();
		
        // Create output variable
		TypedArray<double> outN = solveSoftCrossFields(A_matrix, D_matrix, b_matrix, a_matrix, p, nEdges, nFaces, n, c_matrix, nTotalFaces);
		
        // Set output.
        outputs[0] = outN;
    }
    
	void displayStringOnMATLAB(std::string toprint) {
		return; // Keep quiet for batch testing.
        
        std::ostringstream stream2; stream2 << toprint;
        // Pass stream content to MATLAB fprintf function
        matlabPtr->feval(u"fprintf", 0,
            std::vector<Array>({ factory.createScalar(stream2.str()) }));
        // Clear stream buffer
        stream2.str("");
    }
	
	void displayStreamOnMATLAB(std::ostringstream& stream) {
        return; // Keep quiet for batch testing.
        
        // Pass stream content to MATLAB fprintf function
        matlabPtr->feval(u"fprintf", 0,
            std::vector<Array>({ factory.createScalar(stream.str()) }));
        // Clear stream buffer
        stream.str("");
    }

	void displayCharArrayOnMATLAB(char * buff) {
        return; // Keep quiet for batch testing.
        
		// Pass stream content to MATLAB fprintf function
		std::string somestring(buff);
		std::ostringstream stream2; stream2 << somestring;
		matlabPtr->feval(u"fprintf", 0,
			std::vector<Array>({ factory.createScalar(stream2.str()) }));
		// Clear stream buffer
		stream2.str("");
	}
};