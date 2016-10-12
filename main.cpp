#include "global.hpp"

class Parameters {

	public:
		//declare input
		double alpha;
		double delta_t;		
		double delta_x;	
		double error_max;
		double temperature_initial;
		double temperature_final;
		int N_iter_max;
		int N_nodes;

};

class TDMA_SOLVER_DIRICHLET: public Parameters {

	private:
		const int N_nodes_const = N_nodes;
		std::vector<std::vector<double>> temperature_time;	
		std::vector<double> temperature_spatial;
		double error_computation;

		//compute A constant
		double A_constant_computation (double alpha, double  delta_t, double delta_x) {
			return ( (alpha*delta_t)/pow(delta_x,2) );
		}

		//compute B constant
		double B_constant_computation (double alpha, double  delta_t, double delta_x) {
			return ( 1 + (alpha*delta_t)/pow(delta_x,2) );
		}



	public:
		std::vector<double> temperature_computation (double alpha, 
				double delta_t, 
				double delta_x, 
				double error_max,
				double temperature_initial,
				double temperature_final,
				double N_iter_max,
				double N_nodes) {

			double B = B_constant_computation(alpha, delta_t, delta_x);
			double A = A_constant_computation(alpha, delta_t, delta_x);

			//initialize the vector by initial conditions 
			temperature_spatial.push_back(temperature_initial);
			for (int i=1; i<N_nodes_const-1; i++) {
				temperature_spatial.push_back(0);
			}
			temperature_spatial.push_back(temperature_final);

			//initialize the temperature_time
			temperature_time.push_back(temperature_spatial);

			//looping process
			double time_count = 0;
			while (error_computation > error_max) {


				//clear the temperature_spatial
				temperature_spatial.clear();

				//invoke boundary condition at initial
				temperature_spatial.push_back(temperature_initial);
				
				//compute the 1-st node (initial = 0-th node!)
				double temp_temperature;
				temp_temperature = (temperature_initial - B*temperature_initial)/(-A);
				temperature_spatial.push_back(temp_temperature);
				
				//looping process
				for (int i = 2; i<N_nodes_const; i++) {
					temp_temperature = (temperature_time[time_count][2] + A*temperature_spatial[i-1] - B*temperature_spatial[i])/(-A);
				}




				//invoke boundary condition at final
				temperature_spatial.push_back(temperature_final);

				//add it to temperature_time
				temperature_time.push_back(temperature_spatial);
					
				time_count = time_count + 1;
			};


			std::vector<double> temperature = temperature_spatial;
			return temperature;
		}


};


int main() {

	Parameters par_1;
	//define input of par_1 object
	par_1.alpha			=	0.1;
	par_1.delta_t			=	0.05;
	par_1.delta_x			=	0.1;
	par_1.error_max			=	0.1;
	par_1.temperature_initial	=	0;
	par_1.temperature_final		=	150;
	par_1.N_iter_max		=	100;
	par_1.N_nodes			=	10;

	TDMA_SOLVER_DIRICHLET tdma_1;

	std::vector<double> first = tdma_1.temperature_computation(par_1.alpha, par_1.delta_t, par_1.delta_x, par_1.error_max, par_1.temperature_initial, par_1.temperature_final, par_1.N_iter_max, par_1.N_nodes);

	for (int i = 0; i<par_1.N_nodes; i++) {
		std::cout << first[i] << std::endl;
	
	}





	

}
