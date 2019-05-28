#include <string.h>

#include <mex.hpp>
#include <mexAdapter.hpp>
#include "model.hpp"

//#include <Eigen/Core>
//using namespace matlab::data;
//matlab::data conflict Eigen !!!!!!!
using matlab::mex::ArgumentList;

void show(const model & m)
{
    using namespace  std;
    cout<<"lambda is "<< m.lambda<<endl;
    cout<<"size of mu is "<<m.mu.size()<<endl;
    cout<<"size of nu is "<<m.nu.size()<<endl;
    cout<<"size of pos_mu is "<<m.pos_mu.size()<<endl;
    cout<<"size of pos_nu is "<<m.pos_nu.size()<<endl;
    cout<<"size of M0 is "<<m.M0.size()<<endl;
    double g ;
    m.get_option("CostType",g);
    cout<<"CostType is "<<g<<endl;
    double thr=-1;
    if (m.get_option("thr",thr))
    {
    cout<<"thr is "<<thr<<endl;
    }
    double it;
    m.get_option("min_iter",it);
    cout<<"minimum Sinkhorn iteration is "<<(int)(it)<<endl;

    m.get_option("max_iter",it);
    cout<<"maximum Sinkhorn iteration is "<<(int)(it)<<endl;
    //cout<<"gamma is "<<((m.pskh))<<endl;
    //cout<<"mu=\n" << m.mu<<endl;
    //cout<<"nu=\n" << m.nu<<endl;
    //cout<<"M0=\n" << m.M0<<endl;
}
template<class MATLABARRAY,class EIGENVEC>
void COPYTO(const MATLABARRAY & src, EIGENVEC & v)
{
    std::size_t s = src.getNumberOfElements();
    v.resize(s);
    for(std::size_t i=0;i<s;++i)
    {
	v(i) = src[i];
    }
}
class MexFunction: public matlab::mex::Function{
    
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

    // Factory to create MATLAB data arrays
    matlab::data::ArrayFactory factory;
    // Create an output stream
    std::ostringstream stream;
    int n;
    bool use_log;
    int min_iter;
    int max_iter;
    double cost_type;
    double thr;
    public:
	void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
	    use_log= true;
	    min_iter=20;
	    max_iter=1000;
	    cost_type=1.0;
	    thr = -1;
	    checkArguments(outputs,inputs);
	    using matlab::data::TypedArray;
	    //std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

            //matlabPtr->feval(u"fprintf", 0, 
                //std::vector<matlab::data::Array>({ factory.createScalar("ready") }));
	//stream <<"ready";
	//displayOnMATLAB(stream);
	std::cout<<"ready"<<std::endl;


	    const TypedArray<double> curv = inputs[0];
	    const TypedArray<double> pos = inputs[1];
	    const TypedArray<double> x_in = inputs[2];
	    const TypedArray<double> lambda = inputs[3];
	    n = curv.getNumberOfElements() ;
	    int n2 = pos.getNumberOfElements() ;


	    if(n!=n2)
	    {
            matlabPtr->feval(u"error", 0, 
                std::vector<matlab::data::Array>({ factory.createScalar("the first two inputs must have the same size!") }));
	    }
	    auto dimx= x_in.getDimensions();
	    if(dimx[1]!=4)
	    {
            matlabPtr->feval(u"error", 0, 
                std::vector<matlab::data::Array>({ factory.createScalar("the third input is supposed to be a n * 4 matrix!") }));
	    }
	    int nl = lambda.getNumberOfElements();
	    if(nl!=2)
	    {
            matlabPtr->feval(u"error", 0, 
                std::vector<matlab::data::Array>({ factory.createScalar("the fourth input (lambda) is supposed to contain 2 elements!") }));
	    }
	    int N = dimx[0];
	    typedef Eigen::VectorXd Vec;
	    TypedArray<double> x = factory.createArray<double>(dimx);
	    matlab::data::ArrayDimensions Edim(2);
	    Edim[0]=N; Edim[1] = 1;
	    TypedArray<double> E = factory.createArray<double>(Edim);
	    Vec x_init;
	    x_init.resize(4);
	    Vec x_finish;
	    x_finish.resize(4);
	    //const double & c0 =(curv[0][0]);
	    //const double & p0 =pos[0][0];
	    //const double & l0 =lambda[0][0];
	    Vec c;
	    //c=Eigen::Map<Vec>(&const_cast<double &>(c0),n);
	    Vec p;
	    //p=Eigen::Map<Vec>(&const_cast<double &>(p0),n);
	    Vec l;
	    //l=Eigen::Map<Vec>(&const_cast<double &>(l0),3);
	    COPYTO(curv,c);
	    COPYTO(pos,p);
	    COPYTO(lambda,l);
	    c/=c.sum();
	    model mm(c,p,l,false);
	    mm.print_log = use_log;
	    mm.set_option("min_iter",min_iter);
	    mm.set_option("max_iter",max_iter);
	    if(thr > 0)
	    {
		mm.set_option("thr",thr);
	    }
	    mm.set_option("CostType",cost_type);
	    mm.init_sinkhorn();
	    show(mm);

            //matlabPtr->feval(u"fprintf", 0, 
                //std::vector<matlab::data::Array>({ factory.createScalar("ready") }));
	    for (int term=0;term<N;++term)
	    {
		for(int i=0;i<4;++i)
		{
		    x_init(i) = x_in[term][i];
		}
		E[term] = mm.solve(x_init,x_finish);
		for(int i=0;i<4;++i)
		{
		    x[term][i] = x_finish(i);
		}
	    };

	    outputs[0] = x;
	    if(outputs.size()>=1)
	    {
		outputs[1]=E;
	    };

	};
	void	checkArguments(ArgumentList outputs,ArgumentList inputs)
	{
        //std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        //matlab::data::ArrayFactory factory;
	if (inputs.size()<4)
	{
            matlabPtr->feval(u"error", 0, 
                std::vector<matlab::data::Array>({ factory.createScalar("input: curvature,position,x,lambda,[cost type],[min iter],[max iter],[show_log]") }));
	}

	    using matlab::data::TypedArray;
	    if (inputs.size()>=5)
	    {
		const TypedArray<double> tmp = inputs[4];
		if(! tmp.isEmpty())
		{
		    cost_type = tmp[0];
		    if(cost_type<=0)
		    {
            matlabPtr->feval(u"error", 0, 
                std::vector<matlab::data::Array>({ factory.createScalar("cost_type must be positive!") }));
		    }
		    if(tmp.getNumberOfElements()>=2)
		    {
			thr = tmp[1];
		    }
		}
	    }

	    if (inputs.size()>=6)
	    {
		const TypedArray<double> tmp = inputs[5];
		if(! tmp.isEmpty())
		{
		    min_iter = (int)(tmp[0]);
		    if(min_iter<=0)
		    {
			min_iter = 0;
		    }
		}
	    }

	    if (inputs.size()>=7)
	    {
		const TypedArray<double> tmp = inputs[6];
		if(! tmp.isEmpty())
		{
		    max_iter = (int)(tmp[0]);
		    if(max_iter<=0)
		    {
			max_iter = 0;
		    }
		}
	    }
	    if (inputs.size()>=8)
	    {
		const TypedArray<bool> array_use_log = inputs[7];
		use_log = array_use_log[0];
	    }
        //if (inputs[0].getType() != ArrayType::DOUBLE ||
            //inputs[0].getType() == ArrayType::COMPLEX_DOUBLE)
        //{
            //matlabPtr->feval(u"error", 0, 
                //std::vector<Array>({ factory.createScalar("Input must double array") }));
        //}

        if (outputs.size() >2) {
            matlabPtr->feval(u"error", 0, 
                std::vector<matlab::data::Array>({ factory.createScalar("at most 2 outputs!") }));
        }
	};
    void displayOnMATLAB(std::ostringstream& stream) {
        // Pass stream content to MATLAB fprintf function
        matlabPtr->feval(u"fprintf", 0,
            std::vector<matlab::data::Array>({ factory.createScalar(stream.str()) }));
        // Clear stream buffer
        stream.str("");
    }
};
