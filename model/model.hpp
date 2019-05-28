#pragma once
#include <Eigen/Core>
#include <memory>
#include <map>
#include <string>


//class Sinkhorn_log;
#define SINK Sinkhorn
//#define SINK Sinkhorn
class SINK;
typedef SINK Sinkhorn_solver;

class model 
{
    typedef std::string string;
    typedef std::shared_ptr<Sinkhorn_solver> skh_ptr;
    public:
    typedef std::map<std::string,double> Option;
    typedef std::map<std::string,double> Status;
    typedef Eigen::VectorXd Vec;
    typedef Eigen::MatrixXd Mat;
    typedef Eigen::Ref<Vec> VecRef;
    typedef Eigen::Ref<Mat> MatRef;
    typedef Eigen::Ref<const Vec> VecCRef;
    typedef Eigen::Ref<const Mat> MatCRef;
	model(const Vec& curvature,const Vec& pos,const Vec & lambda,bool init = true);
	void init_sinkhorn();
	static double E_opposite(const Vec & x, Vec & grad);
	static double E_opposite_ad(const Vec & x, Vec & grad);
	//static double E_min_opposite_ad(const Vec & x, Vec & grad);
	double E_Sinkhorn(const Vec & x,Vec & grad);
	double E_Sinkhorn_no_ad(const Vec & x,Vec & grad);
	//static double E_length(const Vec & x, Vec & grad);
	//static double E_length_ad(const Vec & x, Vec & grad);
	static double E_pair_ad(const Vec & x, Vec & grad);
	//static double E_avradj_ad(const Vec & x, Vec & grad);
	//double E_transp(const Vec & x, const Mat & plan, Vec & grad_x, Mat & grad_plan );
	void check_lambda();

	double operator()(const VecCRef & x, VecRef grad);
	//double operator()(const VecCRef & x, const MatCRef & plan, VecRef grad_x, MatRef grad_plan );
	double solve(const Vec& x_init,Vec & x_out);
	double solvek(const Vec& x_init,Vec & x_out);
	void set_option(const string &s,double d);
	bool get_option(const string &s,double & d) const;
	bool get_status(const string &s,double & d) const;

	const Status & status(){return st;};
	Option & options(){return opt;};
	friend void show(const model & m);
	bool print_log;

	skh_ptr pskh;
    private:
	Vec curvature;
	Vec pos;
	Vec lambda;
	Vec mu;
	Vec nu;
	Vec pos_mu;
	Vec pos_nu;
	Mat M0;
	Mat K0;
	Option opt;
	Status st;
};
