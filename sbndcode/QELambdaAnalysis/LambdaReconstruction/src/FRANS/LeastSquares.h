#include <vector>


template <typename T , typename M>
void estimate_coef(std::vector<T> indep_var, std::vector<T> dep_var, M &a1, M &a0, double &r, M &a1err)
{
	int n = indep_var.size();

	M SS_x = 0;
	M SS_y = 0;
  M SS_xy = 0;
  M SS_xx = 0;
  M SS_yy = 0;

	for(int k=0; k<n; k++){
		SS_x+=indep_var[k];
		SS_y+=dep_var[k];
		SS_xy+=indep_var[k]*dep_var[k];
		SS_xx+=indep_var[k]*indep_var[k];
		SS_yy+=dep_var[k]*dep_var[k];
	}

  /*std::cout<< "SS_xy : " << SS_xy <<std::endl;
  std::cout<< "SS_xx : " << SS_xx <<std::endl;
  std::cout<< "SS_x : " << SS_x <<std::endl;
  std::cout<< "SS_y : " << SS_y <<std::endl;*/

	a1 = (n*SS_xy-SS_x*SS_y) / (n*SS_xx-SS_x*SS_x);
 	a0 = (SS_xx*SS_y-SS_x*SS_xy) / (n*SS_xx-SS_x*SS_x);

	//compute errors
	M SS_pred = 0;
	M SS_sigmax = 0;
	M Xmean = SS_x/n;
	for(int k=0; k<n; k++){
		SS_pred+=(dep_var[k]-a1*indep_var[k]-a0)*(dep_var[k]-a1*indep_var[k]-a0);
		SS_sigmax+=(indep_var[k]-Xmean)*(indep_var[k]-Xmean);
	}

	a1err = std::sqrt( SS_pred/( (n-2)*SS_sigmax ) );

	r = (n*SS_xy-SS_x*SS_y)/ std::sqrt( (n*SS_xx-SS_x*SS_x)*(n*SS_yy-SS_y*SS_y) ) ;


}


template <typename X, typename Z>
class Linear_Regression
{
    private:
    X slope;
    X intercept;
		X slope_err;
		double r2;

    public:
    void fit(std::vector<std::vector<Z>> dataset)
    {
        estimate_coef(dataset[0],dataset[1] ,slope, intercept, r2, slope_err);
    }

    Z predict(const Z & test_data)
    {
        return intercept + (slope * test_data);
    }

		X Slope(){
			return slope;
		}

		X SlopeErr(){
			return slope_err;
		}

		X Intercept(){
			return intercept;
		}

		X R2(){
			return r2;
		}


};
