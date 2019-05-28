#include "model.hpp"
#include "Sinkhorn.hpp"
#include <assert.h>
#include <cmath>
#include <vector>
#include <LBFGS.h>
#include <functional>

#include "DEBUG_MESSAGE.hpp"


//#define COST_FUNC(x) ( x )
//#define COST_FUNC(x) ( sqrt(1e-10+x) )
//

    using stan::math::var;
    //using namespace std;
    using Eigen::Dynamic;
    typedef Eigen::Matrix<var,Dynamic,1> aVec;
    typedef Eigen::Matrix<var,Dynamic,Dynamic> aMat;
template<class T>
typename std::function<T(const T&)> get_cost_func(double g)
{
    using namespace std;
    using namespace stan::math;
    if (g == 1)
    {
	return [](const  T & x){ return x;};
    }
    else if (g == 2.0)
    {
	return [] (const T & x){return x*x;};
    }
    else if (g == 0.5) 
    {
	return [] (const T & x){return sqrt(1e-10+x);};
    }
    else if(g>0)
    {
	return [g](const T &x){ return pow(1e-10 +x,g);};
    }
    else
   {
       throw( std::runtime_error("unknown cost function!"));
    }
}


template<class T1,class T2>
typename boost::math::tools::promote_args<T1, T2>::type
smin(const T1 & a,const T2 & b, double k=0.01 )
{
    // a C1 smooth fmin;
    using namespace stan::math;
    using namespace std;
    typedef	typename boost::math::tools::promote_args<T1, T2>::type  T;
    T h = fmax( k-fabs(a-b), 0.0 )/k;
    return fmin( a, b ) - h*h*k*(1.0/4.0);
}
template<class T>
class cost_function_factory
{
    public: 
	const model * m_model_ptr;
	cost_function_factory(const model * m):m_model_ptr(m)
	{
	}
	typename std::function<T(const T&)> get_function()
	{
	    double cost_type=1;
	    m_model_ptr->get_option("CostType",cost_type);
	    auto cost_func = get_cost_func<T>(cost_type);
	    double thr = 1;

	    if( m_model_ptr->get_option("thr",thr))
	    {

    using namespace std;
    using namespace stan::math;
		auto thr_func = [thr](const T& x){return smin(x,thr);};
		return [thr_func,cost_func](const T& x) {return cost_func(thr_func(x));};
	    }
	    return cost_func;
	}
};


template<class V1, class V2>
void subvec_by_ind(const V1 & v1,const V2& index, V1 & out)
{
    std::size_t s = index.size();
    if (out.size()<s) out.resize(s);
    for(std::size_t i=0;i<s;++i)
    {
	out(i) = v1(index[i]);
    }
}


template<class T1>
T1
sabs(const T1 & a, double k=1e-6 )
{
    // a  smoothed fabs;
    using namespace stan::math;
    using namespace std;
    //return 2/k*log(1+exp())
    return a*a/(sqrt(a*a+k));
}



using std::string;
typedef  model::Vec Vec;
model::model(const Vec& curvature,const Vec& pos,const Vec & lambda_in,bool init):curvature(curvature),pos(pos)
{
    //opt["Sinkhorn_eps"] = 1e-3;
    lambda = lambda_in;
    check_lambda();
    opt["CostType"] = 1;
    this->curvature/= curvature.sum();
    if (init)
    {
    init_sinkhorn();
    }
    print_log = false;
};
void model::init_sinkhorn()
{
    // construct mu,nu,M0
    using namespace Eigen;
    std::vector<int> positive;
    std::vector<int> negative;
    static const double srh = 1e-10;
    for (int i =0;i<curvature.size();++i)
    {
	if (curvature(i)>srh)
	    positive.push_back(i);
	else if(curvature(i)<-srh)
	    negative.push_back(i);
    }
    //mu = curvature(positive);//not supported in current version of eigen, only in dev branch
    subvec_by_ind(curvature,positive,mu);
    nu.resize(negative.size()+4);
    //nu.head(negative.size()) = - curvature(negative);
    subvec_by_ind(curvature,negative,nu);
    subvec_by_ind(pos,positive,pos_mu);
    subvec_by_ind(pos,negative,pos_nu);
    nu = -nu;
    for(int i = 0;i<4 ;++i)
    {
	nu(negative.size()+i) = 0.25;
    }
    mu /= mu.sum();
    nu /= nu.sum();
    M0.resize(pos_mu.size(),pos_nu.size()+4);
    M0.setOnes();
    double pm,pn,d;
    for(int i = 0 ;i<pos_mu.size();++i)
    {
	pm = pos_mu(i);
	for(int j = 0;j<pos_nu.size();++j)
	{
	    pn = pos_nu(j);
	    d = std::fmod(10.0+pn-pm,1.0);
	    M0(i,j) = std::min(d,1-d);
	}
    }

    double MM = fmin(M0.maxCoeff(),0.1);
    //pskh = skh_ptr(new Sinkhorn_solver(opt["Sinkhorn_eps"]));
    //gamma = 
    double min_iter=20;
    double max_iter=1000;
    double cost_type=1;
    get_option("min_iter",min_iter);
    get_option("max_iter",max_iter);
    get_option("CostType",cost_type);

    cost_function_factory<double > cfd(this);
    auto COST_FUNC = cfd.get_function();
    //auto COST_FUNC = get_cost_func<double>(cost_type);
    double gamma  =COST_FUNC(MM)/100;
    using std::cout;
    using std::endl;
    cout<<"gamma = " << gamma << endl;
    cout<<"MM = " << MM << endl;

    pskh = skh_ptr(new Sinkhorn_solver(gamma,static_cast<std::size_t>(max_iter),
		static_cast<std::size_t>(min_iter)));
    //if(print_log)
    //std::cout<<"gamma is"<< pskh->gamma<<std::endl;


}

double model::E_opposite_ad(const Vec & x, Vec & grad)
{
    using stan::math::var;
    using namespace stan::math;
    using Eigen::Dynamic;
    typedef Eigen::Matrix<var,Dynamic,1> aVec;
    aVec vx(4);
    grad.resize(4);
    for (int i = 0; i < 4; ++i) {
	vx(i)=x(i);
    }
    aVec l(4);
    //for (int i=0;i<4;++i)
    //{
	//l(i) = fmod(vx((i+1)%4)-vx(i),1.0);
	////grad(i) = -l(i)+l((i-1)%4);
    //}

    for (int i=0;i<4;++i)
    {
	//l(i) = fabs(fmod(vx((i+1)%4)-vx(i),1.0));
	l(i) = (fmod(10.0+vx((i+1)%4)-vx(i),1.0));
	//grad(i) = -l(i)+l((i-1)%4);
    }
	//l(3) = fabs(1 - l(0)-l(1)-l(2));
    var value = 0.5* ((l(0)-l(2))*(l(0)-l(2))+(l(1)-l(3))*(l(1)-l(3)));
    stan::math::set_zero_all_adjoints();
    value.grad();
    for (int i=0;i<4;++i)
    {
	grad(i) =vx(i) .adj();
    }
    recover_memory();
    return value_of(value);
}

double model::E_opposite(const Vec & x, Vec & grad){
    grad.resize(4);
    Vec l(4);
    using namespace std;
    for (int i=0;i<4;++i)
    {
	l(i) = fmod(x((i+1)%4)-x(i),1.0);
	//grad(i) = -l(i)+l((i-1)%4);
    }
    for (int i=0;i<4;++i)
    {
	grad(i) = l((i-1+4)%4)-l((i+1)%4)+l((i+2)%4)-l((i)%4);
    }
    return 0.5* ((l(0)-l(2))*(l(0)-l(2))+(l(1)-l(3))*(l(1)-l(3)));
}
#define SET_UP_M \
\
    aMat M = aMat::Zero(mu.size(),nu.size());\
    for(int i = 0 ;i<pos_mu.size();++i)\
    {\
	pm = pos_mu(i);\
	for(int j = 0;j<pos_nu.size();++j)\
	{\
	    /*M(i,j)=pow(M0(i,j),0.25);*/\
	    M(i,j)=M0(i,j);\
	    /*M(i,j)=sqrt(M0(i,j));*/\
	}\
	for(int j = 0;j<4;++j)\
	{\
	    tmp = sabs(fmod((pm-vx(j)),1.0));\
	    M(i,j+pos_nu.size()) = smin(tmp,1-tmp);\
	    /*M(i,j+pos_nu.size()) = pow(1e-10+M(i,j+pos_nu.size()),0.5)*/;\
	    /*M(i,j+pos_nu.size()) = sqrt(1e-10+M(i,j+pos_nu.size()));*/\
	}\
    };


	    //M(i,j)=(M0(i,j));\
	    //M(i,j+pos_nu.size()) = sqrt(1e-10+M(i,j+pos_nu.size()));

template <typename aVec,typename aMat>
void setup_Matrix(model * self, const Eigen::MatrixXd & M0, const aVec & vx,const Eigen::VectorXd &pos_mu,const Eigen::VectorXd &pos_nu, aMat & M, aMat & K,double g=1)
{

    using namespace stan::math;
    cost_function_factory<double > cfd(self);
    cost_function_factory<var > cfv(self);
    auto COST_FUNC = cfd.get_function();
    auto COST_FUNC_active = cfv.get_function();

    var tmp;
    double pm;
    double gamma = self-> pskh->gamma;
    M.resize(pos_mu.size(),pos_nu.size()+4);\
    K.resize(pos_mu.size(),pos_nu.size()+4);\
    for(int i = 0 ;i<pos_mu.size();++i)\
    {\
	pm = pos_mu(i);\
	for(int j = 0;j<pos_nu.size();++j)\
	{\
	    /*M(i,j)=pow(M0(i,j),0.25);*/\
	    M(i,j)=COST_FUNC(M0(i,j));\
	    K(i,j)= exp(-value_of(M(i,j))/gamma);
	    /*M(i,j)=sqrt(M0(i,j));*/\
	}\
	for(int j = 0;j<4;++j)\
	{\
	    tmp = sabs(fmod((pm-vx(j)),1.0));\
	    M(i,j+pos_nu.size()) = COST_FUNC_active(smin(tmp,1-tmp));\
	    K(i,j+pos_nu.size())= exp(-M(i,j+pos_nu.size())/gamma);
	    /*M(i,j+pos_nu.size()) = pow(1e-10+M(i,j+pos_nu.size()),0.5)*/;\
	    /*M(i,j+pos_nu.size()) = sqrt(1e-10+M(i,j+pos_nu.size()));*/\
	}\
    };

}
double model::E_Sinkhorn_no_ad(const Vec & x,Vec & grad)
{

    double v;
    {
    auto & skh = *pskh;
    //using stan::math::var;
    using namespace std;
    using Eigen::Dynamic;
    typedef Eigen::Matrix<var,Dynamic,1> aVec;
    typedef Eigen::Matrix<var,Dynamic,Dynamic> aMat;
    using namespace::stan::math;
    aVec vx(4);
    grad.resize(4);
    for (int i = 0; i < 4; ++i) {
	vx(i)=x(i);
    }
    var tmp;
    //double pm;
    aMat M;
    aMat K;
    double cost_type = opt["CostType"];

    setup_Matrix(this, M0,vx,pos_mu,pos_nu,M,K,cost_type);
    Eigen::Matrix<double,Dynamic,Dynamic> Mv = value_of(M);
    Eigen::Matrix<double,Dynamic,Dynamic> Kv = value_of(K);


    v = skh(Mv,Kv,mu,nu);
    Eigen::Matrix<double,Dynamic,Dynamic> T = skh.transport_Plan;
    var value = sum(elt_multiply(T,M));
    //var value = sum(M);
    stan::math::set_zero_all_adjoints();
    value.grad();
    for(int i =0;i<4;++i)
    {
	grad(i) = vx(i).adj();
    }
    v = value_of(value);
    skh.clean_init();
#ifndef NDEBUG
    std::cout << "sinkhorn iter: "<<skh.last_iter<<std::endl;
#endif
};
stan::math::recover_memory();


    return v;
}

	    
double model::E_Sinkhorn(const Vec & x,Vec & grad){
    //Sinkhorn_log skh(opt["Sinkhorn_eps"]);
    double v;
    {
    auto & skh = *pskh;
    using stan::math::var;
    using namespace std;
    using Eigen::Dynamic;
    typedef Eigen::Matrix<var,Dynamic,1> aVec;
    typedef Eigen::Matrix<var,Dynamic,Dynamic> aMat;
    using namespace::stan::math;
    aVec vx(4);
    grad.resize(4);
    for (int i = 0; i < 4; ++i) {
	vx(i)=x(i);
    }
    var tmp;
    //double pm;
    aMat M;
    aMat K;
    double cost_type = opt["CostType"];

    setup_Matrix(this, M0,vx,pos_mu,pos_nu,M,K,cost_type);

    var value = skh(M,K,mu,nu);
    //var value = sum(M);
    stan::math::set_zero_all_adjoints();
    value.grad();
    for(int i =0;i<4;++i)
    {
	grad(i) = vx(i).adj();
    }
    v = value_of(value);
    skh.clean_init();
#ifndef NDEBUG
    std::cout << "sinkhorn iter: "<<skh.last_iter<<std::endl;
#endif
};
stan::math::recover_memory();


    return v;
};

double model::E_pair_ad(const Vec & x, Vec & grad)
{
    using stan::math::var;
    using namespace stan::math;
    using Eigen::Dynamic;
    typedef Eigen::Matrix<var,Dynamic,1> aVec;
    aVec vx(4);
    grad.resize(4);
    for (int i = 0; i < 4; ++i) {
	vx(i)=x(i);
    }
    aVec l(4);
    for (int i=0;i<4;++i)
    {
	l(i) = (fmod(10.0+vx((i+1)%4)-vx(i),1.0));
	//l(i) = fabs(fmod(vx((i+1)%4)-vx(i),1.0));
	//grad(i) = -l(i)+l((i-1)%4);
    }
    var value;
    value = 0.5*square(l(0)+l(2)-l(1)-l(3));
    stan::math::set_zero_all_adjoints();
    value.grad();
    for (int i=0;i<4;++i)
    {
	grad(i) =vx(i) .adj();
    }
    recover_memory();
    return value_of(value);
}


void model::check_lambda()
{
    if (lambda.size()!=2)
    {
    std::cout << "the size of lambda should be 2!"<<std::endl;
    }
}
//double model::operator()(const Vec & x, Vec & grad)
double model::operator()(const VecCRef & x, VecRef grad)
{
    assert(lambda.size()==2);
    assert(grad.size()==4);
    grad.setZero();
    if (x(0)>=x(1) || x(1)>=x(2)||x(2)>=x(3) || x(3)>=x(0)+1)
    {

    std::cout<<"!!!!!!!"<<std::endl;
    std::cout<<"x =  "<<x.transpose()<<std::endl;
    std::cout<<"!!!!!!!"<<std::endl;
	return 100;
    };
// v2 : sink op adj pair [avradj]

//#ifndef NDEBUG
    if (print_log)
    {
    std::cout<<"x =  "<<x.transpose()<<std::endl;
    }
//#endif
    Vec g(4);
    double v=0;
    double l;
    l = lambda(0);
    if (l>1e-10)
    {
    double v_op = E_opposite_ad(x,g);
    grad += l*g;
    v += l*v_op;

    if(print_log)
    {
    std::cout<<"v_op =  "<<v_op<<std::endl;
    std::cout<<"grad="<<g.transpose()<<std::endl;
    }
    }



    l = 1.0;
    if (l>1e-10)
    {
    double v_sink = E_Sinkhorn(x,g);
    grad += l*g;
    v+= l*v_sink;


    if(print_log)
    {
    std::cout<<"v_sink =  "<<v_sink<<std::endl;
    std::cout<<"grad="<<g.transpose()<<std::endl;
    std::cout<<"sinkhorn iter = "<<pskh->last_iter<<std::endl;
    }

    }
//#ifndef NDEBUG
    //std::cout<<"sink_iter =  "<<pskh->last_iter<<std::endl;
    //std::cout<<"v_sink =  "<<v_sink<<std::endl;
    //std::cout<<"grad="<<g.transpose()<<std::endl;
//#endif


    l =lambda(1);
    if (l>1e-10)
    {
    double v_pair = E_pair_ad(x,g);

    grad += l*g;
    v+= l*v_pair;
    //DEBUG_MSG("v_l =  "<<v_length);
    //DEBUG_MSG("grad="<<g.transpose());
    
    if(print_log)
    {
    std::cout<<"v_pair =  "<<v_pair<<std::endl;
    std::cout<<"grad="<<g.transpose()<<std::endl;
    }
    }


//#ifndef NDEBUG


    if(print_log)
    {
    std::cout<<"v =  "<<v<<std::endl;
    std::cout<<"|grad|="<<grad.norm()<<std::endl;
    std::cout<<"==========================\n";
    }

//#endif
    return v;
};

template<class T>
struct x_constrants
{
    bool print_log;
    x_constrants(bool b =true):print_log(b){};
    T operator()(const Vec& x,const Vec & dir,T def=1.0)
    {
	//ensure x+s*dir satisfies that
	// x0<x1<x2<x3<x0+1
	double eps = 0.01;
	double step = def;
	for(int i = 0;i<3;++i)
	{
	    if(dir(i)-dir(i+1)>0)
	    step = std::fmin(step,(x(i+1)-x(i)-eps)/(dir(i)-dir(i+1)));
	}
	if(dir(3)-dir(0)>0)
	{
	step = std::fmin(step,(x(0)+1-x(3)-eps)/(dir(3)-dir(0)));
	}
	if(print_log)
	{
	std::cout<<"step= "<<step<<std::endl;
	std::cout<<"x = "<< x.transpose()<<std::endl;

	std::cout<<"dir= "<<dir.transpose()<<std::endl;

	std::cout<<"|dir|= "<<dir.norm()<<std::endl;
	}
	//DEBUG_MSG("step= "<<step);
	//DEBUG_MSG("|dir|= "<<dir.norm());
	
	return T(step);
    }
};




double model::solve(const Vec& x_init,Vec & x){
    assert(x_init.size()==4);
    double fx;
    using namespace LBFGSpp;
    using namespace std;
    LBFGSParam<double> param;
    param.linesearch= LBFGS_LINESEARCH_BACKTRACKING_WOLFE; 
    param.m=50;
    //param.past = 20;
    //param.delta= 1e-4; 

    LBFGSSolver<double> solver(param);
    x_constrants<double> constrants(print_log);
    x = x_init;
    Vec current_x;
    int niter=0;
    Vec grad(4);
    Vec drt;
    double s;
     //gradient method
    
    /*
    for(;niter<param.max_iterations;++niter)
    {
	current_x=x;
	fx = this->operator()(x,grad);
	DEBUG_MSG("begin step");

	std::cout<<"v =  "<<fx<<std::endl;
	std::cout<<"grad="<<grad.transpose()<<std::endl;

	drt = -grad;
	s = constrants(x,drt,std::fmin(1/drt.norm(),1.0));
	DEBUG_MSG("s= "<<s);
	DEBUG_MSG("|drt|= "<<drt.norm());
	//LineSearch<double>::Backtracking(*this,fx,x,grad,s,drt,current_x,param);
	double fx_current=fx;
	for(int i =0;i<10;++i)
	{
	    x = current_x+s*drt;
	fx = this->operator()(x,grad);
	if(fx<fx_current)
	{
	    fx_current = fx;
	    current_x = x;
	}
	else
	{
	    s = 0.5*s;
	}
	}
	DEBUG_MSG("s= "<<s);
	DEBUG_MSG("|drt|= "<<drt.norm());
	DEBUG_MSG("x= "<<x);
	DEBUG_MSG("fx= "<<fx);
	if(drt.norm()<1e-3)
	    break;
    };
    */

    param.max_iterations=100;

    param.pastx = 5;
    param.deltax= 1e-5; 
    niter += solver.minimize(*this,constrants,x,fx);

    //param.max_iterations=150;
    //niter += solver.minimize(*this,constrants,x,fx);
    //niter += solver.minimize(*this,constrants,x,fx);
    //niter += solver.minimize(*this,constrants,x,fx);
    //DEBUG_MSG("restarting.....");
    //niter += solver.minimize(*this,constrants,x,fx);
    //Vec g;
    //int niter =0;
	//while(niter<50)
	//{
	//niter+= solver.minimize(*this,x,fx) ;
	//this->operator()(x,g);
	//if(g.norm()<0.01)
	    //break;
	
	//}
	Vec g(4);
	this->operator()(x,g);
	cout<<"|grad*|= "<<g.norm()<<endl;
    cout<<"niter_LBFGS="<<niter<<endl;
    return fx;
    
};
void model::set_option(const string &s,double d){
    opt[s] = d;
};
bool model::get_option(const string &s,double & d) const {
    auto f=opt.find(s);
    bool b = f != opt.end();
    if (b)
    {
	d = f->second;
    }
    return b;
};
bool model::get_status(const string &s,double & d) const {

    auto f=st.find(s);
    bool b = f!= st.end();
    if (b)
    {
	d = f->second;
    }
    return b;
};

