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
		double error_computation;

		//compute A constant
		double A_constant_computation (double alpha, double  delta_t, double delta_x) {
			return ( (alpha*delta_t)/pow(delta_x,2) );
		}

		//compute B constant
		double B_constant_computation (double alpha, double  delta_t, double delta_x) {
			return ( 1 + (alpha*delta_t)/pow(delta_x,2) );
		}
		
		//computa P and Q constant
		double P_init_comp (double b_const, double c_const) {
			return ( (c_const)/(b_const) );
		}

		double P_comp (double a_const, double b_const, double c_const, double P_before) {
			return ( (c_const)/(b_const - a_const*P_before) );	
		}

		double Q_init_comp (double b_const, double d_const) {
			return ( d_const/b_const );
		}

		double Q_comp (double a_const, double b_const, double d_const, double Q_before, double P_before) {
			return ( (d_const - a_const*Q_before)/(b_const - a_const*P_before) );
		}


	public:
		std::vector<double> temperature_computation (double alpha, 
				double delta_t, 
				double delta_x, 
				double error_max,
				double temperature_initial,
				double temperature_final,
				int N_iter_max,
				int N_nodes_const) {

			std::vector<double> temperature_spatial(N_nodes_const);
			std::vector<std::vector<double>> temperature_time(N_iter_max, std::vector<double>(N_nodes_const));

			double B = B_constant_computation(alpha, delta_t, delta_x);
			double A = A_constant_computation(alpha, delta_t, delta_x);


			//initialize the vector by initial conditions 
			temperature_spatial[0] = temperature_initial;
			for (int i=1; i<N_nodes_const-1; i++) {
				temperature_spatial[i] = 0.;
			}
			temperature_spatial[N_nodes_const-1] = temperature_final;

			//initialize the temperature_time
			std::cout << "step[0]: ";
			for (int i=0; i<N_nodes_const; i++) {
				temperature_time[0][i] = temperature_spatial[i];
				std::cout << temperature_time[0][i] << " ";
			}
			std::cout << std::endl;

			int time_count = 0;
			double error_computation = 0.;
			std::vector<double> P(N_nodes_const);
			std::vector<double> Q(N_nodes_const);

			//looping process
			do {

				//compute P and Q constant
				//first compute P, the initial, after that the whole P for N_nodes_const
				P[0]	=	P_init_comp(B, A);
				for (int i = 1; i<N_nodes_const; i++) {
					P[i]	=	P_comp(A, B, A, P[i-1]);
				}

				//after that, compute Q from the initial until N_nodes_const
				Q[0]	=	Q_init_comp(B, temperature_initial);
				for (int i = 1; i<N_nodes_const; i++) {
					Q[i]	=	Q_comp(A, B, temperature_time[time_count][i], Q[i-1], P[i-1]);
				}
	
				//backward substitution for temperature
				temperature_spatial[N_nodes_const - 1]	=	temperature_final;//Q[N_nodes_const - 1];
				for (int i = N_nodes_const - 2; i >= 0; i-- ) {
					temperature_spatial[i]	=	Q[i] - P[i]*temperature_spatial[i+1];
				}

				time_count = time_count + 1;
				std::cout << "step[" << time_count<< "]: "; 
				for (int i = 0; i<N_nodes_const; i++) {
					std::cout << temperature_spatial[i] << " ";
				}
				std::cout << std::endl;


				for (int i=0; i<N_nodes_const; i++) {
					temperature_time[time_count][i] = temperature_spatial[i];
				}

				//compute error
				for (int i = 0; i <N_nodes_const; i++) {
					double temp = abs(temperature_time[time_count][i] - temperature_time[time_count - 1][i]);
					error_computation = (error_computation + temp);		
				}
				error_computation = error_computation/N_nodes_const;
	
				if (time_count > 2) {
					std::cout << "GAGAL" << std::endl;
					break;
				}


			} while (error_computation > error_max);

		return temperature_spatial;










/*

			//initialize the vector by initial conditions 
			temperature_spatial[0] = temperature_initial;
			for (int i=1; i<N_nodes_const-1; i++) {
				temperature_spatial[i] = 0.;
			}
			temperature_spatial[N_nodes_const-1] = temperature_final;
			//temperature_spatial.push_back(temperature_initial);
			//for (int i=1; i<N_nodes_const-1; i++) {
			//	temperature_spatial.push_back(0);
			//}
			//temperature_spatial.push_back(temperature_final);


			//initialize the temperature_time
			//temperature_time.push_back(temperature_spatial);
			for (int i=0; i<N_nodes_const; i++) {
				temperature_time[0][i] = temperature_spatial[i];
			}

			//looping process
			do {

				//clear the temperature_spatial
				//temperature_spatial.clear();
				for (int i=0; i<N_nodes_const; i++) {
					temperature_spatial[i] = 0.;
				}

				//invoke boundary condition at initial
				//temperature_spatial.push_back(temperature_initial);
				temperature_spatial[0] = temperature_initial;

				//invoke boundary condition at final
				//temperature_spatial.push_back(temperature_final);
				temperature_spatial[N_nodes_const-1] = temperature_final;

				//compute the 1-st node (initial = 0-th node!)
				//temperature_spatial.push_back(temp_temperature);
				temperature_spatial[1] = (temperature_initial - B*temperature_initial)/(-A);

				//looping process
				for (int i = 1; i<N_nodes_const-2; i++) {
					//temperature_spatial.push_back(temp_temperature);
					temperature_spatial[i+1] = (temperature_time[time_count][i] + A*temperature_spatial[i-1] - B*temperature_spatial[i])/(-A);					
				}

				time_count = time_count + 1;
				//add it to temperature_time
				//temperature_time.push_back(temperature_spatial);
				for (int i=0; i<N_nodes_const; i++) {
					temperature_time[time_count][i] = temperature_spatial[i];
					std::cout << temperature_spatial[i] << "\t";
				}
				std::cout << "\n";


				for (int i = 0; i <N_nodes_const; i++) {
					double temp = abs(temperature_time[time_count][i] - temperature_time[time_count - 1][i]);
					error_computation = (error_computation + temp);		
				}
				error_computation = error_computation/N_nodes_const;
			
			
				if (time_count > 10) {
					std::cout << "GAGAL" << std::endl;
					break;
				}

			} while (error_computation > error_max);

			std::vector<double> temperature = temperature_spatial;
			return temperature;
			*/
		}
};


int main() {

	Parameters par_1;
	//define input of par_1 object
	par_1.alpha			=	0.1;
	par_1.delta_t			=	0.05;
	par_1.delta_x			=	0.1;
	par_1.error_max			=	0.1;
	par_1.temperature_initial	=	1;
	par_1.temperature_final		=	150;
	par_1.N_iter_max		=	100;
	par_1.N_nodes			=	10;

	TDMA_SOLVER_DIRICHLET tdma_1;

	std::vector<double> first = tdma_1.temperature_computation(par_1.alpha, par_1.delta_t, par_1.delta_x, par_1.error_max, par_1.temperature_initial, par_1.temperature_final, par_1.N_iter_max, par_1.N_nodes);

	for (int i = 0; i<par_1.N_nodes; i++) {
		std::cout << first[i] << "\t" ;
	
	}
	std::cout << std::endl;


}
