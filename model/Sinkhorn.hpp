#pragma once
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <stan/math.hpp>
#include <cmath>
#include <assert.h>
#include "DEBUG_MESSAGE.hpp"

# define MAX_ITER_DEFAULT 1000
namespace stan{
namespace math{
template <typename Ta, int Ra, int Ca, typename Tb, int Cb>
class multiply_mat_buffer_vari : public vari {
 public:
  int A_rows_;
  int A_cols_;
  int B_cols_;
  int A_size_;
  int B_size_;
  double* Ad_;
  double* Bd_;
  vari** variRefA_;
  vari** variRefB_;
  vari** variRefAB_;

  /**
   * Constructor for multiply_mat_vari.
   *
   * All memory allocated in
   * ChainableStack's stack_alloc arena.
   *
   * It is critical for the efficiency of this object
   * that the constructor create new varis that aren't
   * popped onto the var_stack_, but rather are
   * popped onto the var_nochain_stack_. This is
   * controlled to the second argument to
   * vari's constructor.
   *
   * @param A matrix
   * @param B matrix
   */
  multiply_mat_buffer_vari(const Eigen::Matrix<Ta, Ra, Ca>& A,
                    const Eigen::Matrix<Tb, Ca, Cb>& B, double * vA=NULL,double * vB=NULL)
      : vari(0.0),
        A_rows_(A.rows()),
        A_cols_(A.cols()),
        B_cols_(B.cols()),
        A_size_(A.size()),
        B_size_(B.size()),
        variRefA_(
            ChainableStack::instance().memalloc_.alloc_array<vari*>(A_size_)),
        variRefB_(
            ChainableStack::instance().memalloc_.alloc_array<vari*>(B_size_)),
        variRefAB_(ChainableStack::instance().memalloc_.alloc_array<vari*>(
            A_rows_ * B_cols_)) {
	    using Eigen::Map;
	    using Eigen::MatrixXd;



	    if (NULL == vA)
	    {
		Ad_=ChainableStack::instance().memalloc_.alloc_array<double>(A_size_);
		for (size_type i = 0; i < A.size(); ++i) {
		    variRefA_[i] = A.coeffRef(i).vi_;
		    Ad_[i] = A.coeffRef(i).val();
		}
	    }
	    else
	    {
	    for (size_type i = 0; i < A.size(); ++i) {
		    variRefA_[i] = A.coeffRef(i).vi_;
	    }
	    Ad_ = vA;}



	    if(NULL == vB)
	    {
		Bd_=(ChainableStack::instance().memalloc_.alloc_array<double>(B_size_));
	    for (size_type i = 0; i < B.size(); ++i) {
		variRefB_[i] = B.coeffRef(i).vi_;
		Bd_[i] = B.coeffRef(i).val();
	    }
	    }
	    else
	    {
	    for (size_type i = 0; i < B.size(); ++i) {
		variRefB_[i] = B.coeffRef(i).vi_;
	    }
		Bd_ = vB;}
	    MatrixXd AB = Map<MatrixXd>(Ad_, A_rows_, A_cols_)
		* Map<MatrixXd>(Bd_, A_cols_, B_cols_);
	    for (size_type i = 0; i < AB.size(); ++i)
		variRefAB_[i] = new vari(AB.coeffRef(i), false);
  }

  virtual void chain() {
    using Eigen::Map;
    using Eigen::MatrixXd;
    MatrixXd adjAB(A_rows_, B_cols_);
    MatrixXd adjA(A_rows_, A_cols_);
    MatrixXd adjB(A_cols_, B_cols_);

    for (size_type i = 0; i < adjAB.size(); ++i)
      adjAB(i) = variRefAB_[i]->adj_;
    adjA.noalias() = adjAB * Map<MatrixXd>(Bd_, A_cols_, B_cols_).transpose();
    adjB.noalias() = Map<MatrixXd>(Ad_, A_rows_, A_cols_).transpose() * adjAB;
    for (size_type i = 0; i < A_size_; ++i)
      variRefA_[i]->adj_ += adjA(i);
    for (size_type i = 0; i < B_size_; ++i)
      variRefB_[i]->adj_ += adjB(i);
  }
};


template < int Ra, int Ca ,int Cb>
inline  Eigen::Matrix<double, Ra, Cb> 
multiply_b(const Eigen::Matrix<double, Ra, Ca>& A,
         const Eigen::Matrix<double, Ca, Cb>& B,double *vA=NULL,double * vB=NULL) {
  //check_multiplicable("multiply", "A", A, "B", B);
  //check_not_nan("multiply", "A", A);
  //check_not_nan("multiply", "B", B);

  //// Memory managed with the arena allocator.
  //multiply_mat_buffer_vari<Ta, Ra, Ca, Tb, Cb>* baseVari
      //= new multiply_mat_buffer_vari<Ta, Ra, Ca, Tb, Cb>(A, B,vA,vB);
  //Eigen::Matrix<var, Ra, Cb> AB_v(A.rows(), B.cols());
  //for (size_type i = 0; i < AB_v.size(); ++i) {
    //AB_v.coeffRef(i).vi_ = baseVari->variRefAB_[i];
  //}
  //return AB_v;
  return A*B;
}





template <typename Ta, int Ra, int Ca, typename Tb, int Cb>
inline typename boost::enable_if_c<boost::is_same<Ta, var>::value
                                       || boost::is_same<Tb, var>::value,
                                   Eigen::Matrix<var, Ra, Cb> >::type
multiply_b(const Eigen::Matrix<Ta, Ra, Ca>& A,
         const Eigen::Matrix<Tb, Ca, Cb>& B,double *vA=NULL,double * vB=NULL) {
  check_multiplicable("multiply", "A", A, "B", B);
  check_not_nan("multiply", "A", A);
  check_not_nan("multiply", "B", B);

  // Memory managed with the arena allocator.
  multiply_mat_buffer_vari<Ta, Ra, Ca, Tb, Cb>* baseVari
      = new multiply_mat_buffer_vari<Ta, Ra, Ca, Tb, Cb>(A, B,vA,vB);
  Eigen::Matrix<var, Ra, Cb> AB_v(A.rows(), B.cols());
  for (size_type i = 0; i < AB_v.size(); ++i) {
    AB_v.coeffRef(i).vi_ = baseVari->variRefAB_[i];
  }
  return AB_v;
}


}// namespace math
}// namespace stan




struct Sinkhorn
{
    typedef Eigen::VectorXd Vec;
    Sinkhorn(double gamma=1e-2,std::size_t it = MAX_ITER_DEFAULT,std::size_t min_iter=10):max_iter(it),min_iter(min_iter),gamma(gamma){clean_init(); srh_err = 1e-3;};
    //template <class T>
	//T operator()(Eigen::SparseMatrix<T,Eigen::Dynamic,Eigen::Dynamic> & K,Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & M,const Vec & mu, const Vec & nu)
	//{};
    template <class T>
	T operator()(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & M,const Vec & mu, const Vec & nu)
	{
	    using namespace stan::math;
	    using namespace std;
	    using namespace Eigen;
	    typedef typename Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> TMat;
	    TMat K1 = exp(-M/gamma);
	    return  this->operator()(M,K1,mu,nu);
	};
    template <class T>
	T operator()(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & M,Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & K1,const Vec & mu_in, const Vec & nu_in)
	{
	    typedef typename Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> TMat;
	    typedef typename Eigen::Matrix<T,Eigen::Dynamic,1> TVec;
	    using namespace stan::math;
	    using namespace std;
	    using namespace Eigen;
	    assert(M.rows()==mu_in.size());
	    assert(M.cols()==nu_in.size());
	    if(u0.size()!=mu_in.size())
	    {
		u0 = Vec::Ones(mu_in.size())/mu_in.size();
	    }
	    Mv =stan::math::value_of(M);
	    double tolarence = fmin(gamma/(8*((Mv.cwiseAbs().maxCoeff()))),5e-3);
	    double total_mass = sum(mu_in);
	    Vec mu = mu_in/total_mass;
	    Vec nu = nu_in/total_mass;
	    //nu = nu/total_mass;
	    //T min_M = min(M);
	    //double min_M=0;
	    //cout << "gamma"<<gamma<<endl;
	    //cout<< "mu"<<mu<<endl;
	    double tau = -0.5;
	    //Mat mm = Mat::Constant(M.rows(),M.cols(),min_M);
	    //TMat K1 = exp(-(M.array()-min_M).matrix()/gamma);
	    //TMat K1 = exp(-M/gamma);
	    TMat TK1 = transpose(K1);
	    K = value_of(K1);
	    Kt = value_of(TK1);
	    TVec u(u0.size());
	    for (int i = 0 ;i< u0.size();++i)
	    {
		u(i) = u0(i);
	    }
	    TVec v ;
	    TVec Ktu;
	    TVec Kv;
	    for(last_iter = 0;last_iter<max_iter;++last_iter)
	    {
		//Ktu=stan::math::multiply(transpose(K1),u);
		Ktu=stan::math::multiply_b(TK1,u,Kt.data());
		//v = v*tau + (1-tau)*stan::math::elt_divide(nu,Ktu);
		//v = stan::math::elt_divide(nu,add(Ktu,1e-16));
		v = stan::math::elt_divide(nu,Ktu);
		//Kv= stan::math::multiply(K1,v);
		Kv= stan::math::multiply_b(K1,v,K.data());
		//u = u*tau + (1-tau)*stan::math::elt_divide(mu,Kv);
		u = stan::math::elt_divide(mu,Kv);
		//u = stan::math::elt_divide(mu,add(Kv,1e-16));
		if (last_iter>=min_iter 
			&&  (elt_multiply(value_of(v),value_of(Ktu))-nu).norm()<tolarence 
			&& (elt_multiply(value_of(u),value_of(Kv))-mu).norm()<tolarence
			)
		{
		    // checking one side is enough
		    break;
		}
	    }
	    u0 = value_of(u);
	    //cout<<"K1= "<<K1<<endl;
	    //cout<<u<<endl;
	    //cout<<v<<endl;
	    TMat transport_plan =diag_post_multiply(diag_pre_multiply(u,K1),v);
	    this->transport_Plan = stan::math::value_of(transport_plan);
	    //cout<< transport_plan<<endl;
	    //cout<<entropy(transport_plan)<<endl;
	    //
	    //cout<<"transport plan:"<<endl;
	    //cout<<transport_plan<<endl;
	    //T result =sum( elt_multiply(value_of(transport_plan),M)) - gamma*value_of(entropy(transport_plan));
	    //T result = dot_product((u),multiply(elt_multiply(value_of(K1),M),v)) - gamma*value_of(entropy(transport_plan));
	    T result = (dot_product(u,multiply(elt_multiply(K1,M),v)));
	    //T result =(dot_product(transport_plan.data(),M.data(),M.size()));
		    //- total_mass*gamma*(entropy(transport_plan)));
	    //T result = sum(elt_multiply(transport_plan,M)) - gamma*(entropy(transport_plan));
	    //cout<<"max plan"<<value_of(transport_plan).maxCoeff()<<endl;
	    //cout<<"max M"<<Mv.maxCoeff()<<endl;
	    //cout<<"max u"<<endl;
	    //cout<<value_of(u).maxCoeff()<<endl;
	    //cout<<"max v"<<endl;
	    //cout<<value_of(v).maxCoeff()<<endl;
	    //cout<<"entropy="<<entropy(transport_plan)<<endl;
	    //cout<<"v1 = "<<sum(elt_multiply(value_of(transport_plan),Mv))<<endl;
	    //cout<<"v2 = "<<gamma*entropy(transport_plan)<<endl;
	    PrimalE = value_of(result);
	    
	    Vec alpha = log(value_of(u));
	    Vec beta = log(value_of(v));
	    DualE =( gamma*(dot_product(alpha,mu)+dot_product(nu,beta)));
	    //DualE =total_mass*( gamma*(dot_product(alpha,mu)+dot_product(nu,beta))
	    //- gamma*(sum(value_of(transport_plan))));
	    //cout<<"err = "<<(elt_multiply(value_of(v),value_of(Ktu))-nu).norm()<<endl;
	    //
#ifndef NDEBUG
	    //cout << "|mu|,|nu| " << sum(mu)<<" , "<<sum(nu)<<endl;
	    //cout<<"sum of plan "<< sum(value_of(transport_plan))<<endl;
		cout<<"primal = \t"<<PrimalE<<endl;
		cout<<"dual   = \t"<<DualE<<endl;
#endif

	    return result;
	};
    template <class T>
	T entropy(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & m)
	{
	    using namespace stan::math;
	    using namespace std;
	    using namespace Eigen;
	    return - (m.array() * ((m.array()+1e-20).log()-1)).sum();
	}
    void clean_init(){u0.resize(0);};
    std::size_t max_iter;
    std::size_t min_iter;
    std::size_t last_iter;
    double gamma;
    Vec u0;
    double srh_err;
    double PrimalE;
    double DualE;
    Eigen::MatrixXd Mv;
    Eigen::MatrixXd K;
    Eigen::MatrixXd Kt;
    Eigen::MatrixXd transport_Plan;
};

struct Sinkhorn_log:public Sinkhorn
{
    // this one is not fast enough
    Sinkhorn_log(double gamma=1e-2,std::size_t it = MAX_ITER_DEFAULT,std::size_t min_iter=10):Sinkhorn(gamma,it,min_iter){};
    template <class T1,class T2,class T3>
	typename boost::math::tools::promote_args<T1, T2,T3>::type operator()(const Eigen::Matrix<T1,Eigen::Dynamic,Eigen::Dynamic> & M,const Eigen::Matrix<T2,Eigen::Dynamic,1> & mu, const Eigen::Matrix<T3,Eigen::Dynamic,1>& nu){

	    // the actual scaling factor u_t = u.*exp(alpha/gamma)
	    // v_t = v.*exp(beta/gamma)
	    typedef typename boost::math::tools::promote_args<T1, T2,T3>::type T;
	    typedef typename Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> TMat;
	    typedef typename Eigen::Matrix<T,Eigen::Dynamic,1> TVec;
	    typedef typename Eigen::VectorXd Vec;
	    typedef typename Eigen::MatrixXd Mat;
	    using namespace stan::math;
	    using namespace std;
	    using namespace Eigen;
	    assert(M.rows()==mu.size());
	    assert(M.cols()==nu.size());
	    //assert(alpha.size()==mu.size());
	    //assert(beta.size()==nu.size());
	    Mat Mv =stan::math::value_of(M);
	    //if(alpha0.size()!=mu.size()) { alpha0 = Mv.rowwise().minCoeff(); };
	    //Mat Ma = Mv.colwise()-alpha0;
	    //if(beta0.size()!=nu.size()) { beta0 = Ma.colwise().minCoeff(); };
	    if(alpha0.size()!=mu.size()) { alpha0 = Vec::Zero(mu.size()); };
	    if(beta0.size()!=nu.size()) { beta0 = Vec::Zero(nu.size());};
	    double tolarence = gamma/(8*((Mv.cwiseAbs().maxCoeff())));
	    //double tolarence = 1e-5;
	    double Tau=1e10;
	    double tau = -0.05;

	    TVec alpha(mu.size());
	    TVec beta(nu.size());
	    for(int i = 0;i<mu.size();++i){alpha(i)=alpha0(i);}
	    for(int i = 0;i<nu.size();++i){beta(i)=beta0(i);}
	    //copy vectors... there should be better ways.

	    TVec u= TVec::Ones(mu.size());
	    TVec v= TVec::Ones(nu.size());
	    TMat K = getK(M,value_of(alpha),value_of(beta));


		    #ifndef NDEBUG
		//cout<<"K=\n"<<K<<endl;
#endif

	    TVec Ktu;
	    TVec Kv;
	    bool cirt=false;
	    Mat transport_plan(M.rows(),M.cols());
	    Vec valpha;
	    Vec vbeta;
	    for(last_iter = 0;last_iter<max_iter;++last_iter)
	    {
#ifndef NDEBUG
		//cout<<"iter="<<last_iter<<endl;
		//cout<<"u=\n"<<u<<endl;
#endif 
		Ktu=multiply(transpose(K),u);
		v = elt_divide(nu,add(Ktu,1e-16));
		//v =tau*v+(1-tau)* elt_divide(nu,add(Ktu,1e-16));
		Kv= multiply(K,v);
		//u =tau*u+(1-tau)* elt_divide(mu,add(Kv,1e-16));
		u =elt_divide(mu,add(Kv,1e-16));
		if(value_of(u).cwiseAbs().maxCoeff()>Tau ||
			value_of(v).cwiseAbs().maxCoeff()>Tau)
		{
		    alpha += gamma*(log(u));
		    beta += gamma *(log(v));
		    u= TVec::Ones(mu.size());
		    v= TVec::Ones(nu.size());
		    K = getK(M,value_of(alpha),value_of(beta)); // also the current transport plan
		    DEBUG_MSG("absorb\n");
		 if(last_iter%10==0)
		 {
		    transport_plan = value_of(K);
		 };

		}
		else
		{

		//transport_plan =
		    //diag_post_multiply(diag_pre_multiply(value_of(u),value_of(K)),value_of(v));
		 if(last_iter%10==0)
		 {
		  transport_plan = getK(Mv,value_of(alpha),value_of(beta),value_of(u),value_of(v));
		 };
		}

#ifndef NDEBUG
		    //cout<<"transport_plan=\n"<<transport_plan<<endl;
#endif

		//cout<<"transport plan:"<<endl;
		//cout<<transport_plan<<endl;
		//cout << transport_plan.rowwise().sum()<<endl;
		//cout << "mu = "<<mu<<endl;
		//cout << "err=" <<(transport_plan.rowwise().sum()-mu).norm()<<endl;

		//cout<<"======="<<endl;
		//cout << transport_plan.colwise().sum()<<endl;
		//cout << "nu = "<<nu<<endl;

		//cout << "err=" <<
		//(transport_plan.colwise().sum().transpose()-nu).norm()
		//<<endl;

		//cout<<"===================="<<endl;
		if(last_iter%10==0)
		{
		valpha = value_of(alpha)+gamma*(log(value_of(u)));
		vbeta  = value_of(beta)+gamma*(log(value_of(v)));

		PrimalE = value_of(sum(elt_multiply(Mv,transport_plan))-gamma*entropy(transport_plan));
		DualE = dot_product(valpha,mu)+dot_product(nu,vbeta)-gamma*sum(transport_plan);
#ifndef NDEBUG
		cout<<"primal = \t"<<PrimalE<<endl;
		cout<<"dual = \t"<<DualE<<endl;
#endif

		cirt = (last_iter >= min_iter) &&(
		    (fabs(PrimalE-DualE)< srh_err) ||
		    ( ((transport_plan.rowwise().sum()-mu).norm()<tolarence) &&
		    ((transport_plan.colwise().sum().transpose()-nu).norm()<tolarence)));
		if(cirt)
		{
		    break;
		}
		};
	    }
	    alpha += gamma*(log(u));
	    beta += gamma *(log(v));
	    //u= TVec::Ones(mu.size());
	    //v= TVec::Ones(nu.size());
	    K = getK(M,alpha,beta); 
	    alpha0 = value_of(alpha);
	    beta0 = value_of(beta);
	    T result = sum(elt_multiply(M,K))-gamma*entropy(K);
	    PrimalE = value_of(result);
	    DualE = (dot_product(value_of(alpha),mu)+dot_product(nu,value_of(beta)))
	    - gamma*sum(transport_plan) ;
	    //cout<<"sinkhorn iter = "<<last_iter<<endl;

#ifndef NDEBUG
		cout<<"primal = \t"<<PrimalE<<endl;
		cout<<"dual = \t"<<DualE<<endl;
		//cout<<"err1 = "<<(transport_plan.rowwise().sum()-mu).norm()<<endl;
		//cout<<"err2 vec = "<<transport_plan.colwise().sum().transpose()-nu<<endl;
		//cout<<"err2 = "<<(transport_plan.colwise().sum().transpose()-nu).norm()
//<<endl;
#endif
	    return result;

	}
    ;
    template<class T1,class T2,class T3,int R,int C>
	Eigen::Matrix<typename boost::math::tools::promote_args<T1, T2,T3>::type, R, C> getK(const Eigen::Matrix<T3, R, C> & M, const Eigen::Matrix<T1, R, 1> &alpha, const Eigen::Matrix<T2, C, 1> &beta )
	{
	    typedef typename boost::math::tools::promote_args<T1, T2,T3>::type T;
	    typedef typename Eigen::Matrix<T,R,C> TMat;
	    TMat M1 = M.colwise() -alpha;
	    TMat M2 = M1.rowwise() - beta.transpose();

#ifndef NDEBUG
	    //std::cout<<"M2\n"<<M2<<std::endl;
#endif
	    return exp(-1/gamma*(M2));
	}

    template<class T1,class T2,class T3,class T4,class T5,int R,int C>
	Eigen::Matrix<typename boost::math::tools::promote_args<T1, T2,T3,T4,T5>::type, R, C> getK(const Eigen::Matrix<T3, R, C> & M, const Eigen::Matrix<T1, R, 1> &alpha, const Eigen::Matrix<T2, C, 1> &beta, const Eigen::Matrix<T4, R, 1> & u, const Eigen::Matrix<T5, C, 1> &v)
	{
	    typedef typename boost::math::tools::promote_args<T1, T2,T3,T4,T5>::type T;
	    typedef typename Eigen::Matrix<T,R,C> TMat;
	    using namespace stan::math;
	    //TMat M1 = M.colwise() -alpha;
	    //TMat M2 = M1.rowwise() - beta.transpose();
	    TMat M1 = ((-1/gamma*((M.colwise() -alpha).rowwise()- beta.transpose())).colwise()+log(u) ).rowwise()+log(v).transpose();


#ifndef NDEBUG
	    //std::cout<<"M2\n"<<M2<<std::endl;
#endif
	    return stan::math::exp(M1);
	}

    //std::size_t max_iter;
    //double gamma;
    //Vec u0;
    //Vec v0;
    Vec alpha0;
    Vec beta0;
    void clean_init()
    {
	Sinkhorn::clean_init();
	alpha0.resize(0);
	beta0.resize(0);
    }
    void set_init(const Vec & alpha,const Vec & beta)
    {
	alpha0=alpha;
	beta0=beta;
    }


};

