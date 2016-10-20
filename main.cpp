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

		//compute B constant WHERE AS B = 1 + 2A
		double B_constant_computation (double alpha, double  delta_t, double delta_x) {
			return ( 1 + (2*alpha*delta_t)/pow(delta_x,2) );
		}
		
		//computa P and Q constant
//		double P_init_comp (double b_const, double c_const) {
//			return ( c_const/b_const );
//		}

//		double P_comp (double a_const, double b_const, double c_const, double P_before) {
//			return ( (c_const)/(b_const - a_const*P_before) );	
//		}

//		double Q_init_comp (double b_const, double d_const) {
//			return ( d_const/b_const );
//		}

//		double Q_comp (double a_const, double b_const, double d_const, double Q_before, double P_before) {
//			return ( (d_const - a_const*Q_before)/(b_const - a_const*P_before) );
//		}


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

			//make the initial matrix
			std::vector<std::vector<double>> the_matrix(N_nodes_const,std::vector<double>(N_nodes_const));

			for (int j=0; j<N_nodes_const; j++) {	//vertical -> rows
				for (int i=0; i<N_nodes_const; i++) {	//horizontal -> columns

					if (j == 0) {
						if (i == 0) {
							the_matrix[j][i]	=	1;
						}
						else {
							the_matrix[j][i]	=	0;
						}
					}

					else if (j == N_nodes_const - 1) {
						if (i == N_nodes_const - 1) {
							the_matrix[j][i]	=	1;
						}
						else {
							the_matrix[j][i]	=	0;
						}
					}
					
					else {
						if ( i == j ) {
							the_matrix[j][i]	=	B;
						}
						else if (i == j-1 || i == j+1) {
							the_matrix[j][i]	=	-A;
						}
						else {
							the_matrix[j][i] 	=	0;
						}
					}


				}
			}
/*			for (int j=0; j<N_nodes_const; j++) {
				for (int i=0; i<N_nodes_const; i++) {
					std::cout << the_matrix[j][i] << "\t";
				}
				std::cout << std::endl;
			}
*/

			int time_count = 0;
			double error_computation = 0.;
			//std::vector<double> P(N_nodes_const);
			//std::vector<double> Q(N_nodes_const);

			//looping process
			do {

				//compute P and Q constant
				//first compute P, the initial, after that the whole P for N_nodes_const
				//P[0]	=	P_init_comp(B, -A);
				//for (int i = 1; i<N_nodes_const; i++) {
				//	P[i]	=	P_comp(-A, B, -A, P[i-1]);
				//}

				//after that, compute Q from the initial until N_nodes_const
				//Q[0]	=	Q_init_comp(B, temperature_initial);
				//for (int i = 1; i<N_nodes_const; i++) {
				//	Q[i]	=	Q_comp(A, B, temperature_time[time_count][i], Q[i-1], P[i-1]);
				//}

				//backward substitution for temperature
				//temperature_spatial[N_nodes_const - 1]	=	temperature_final;//Q[N_nodes_const - 1];
				//for (int i = N_nodes_const - 2; i >= 0; i-- ) {
				//	temperature_spatial[i]	=	Q[i] - P[i]*temperature_spatial[i+1];
				//}


				// ------------------------ do the gauss method ------------------
				// first, make a bottom triangle matrix
				for (int j = 1; j<N_nodes_const-1; j++) {
					double temp_const = the_matrix[j][j-1]/the_matrix[j-1][j-1];
					for (int i=0; i<N_nodes_const; i++) {
						the_matrix[j][i] = the_matrix[j][i] - the_matrix[j-1][i]*temp_const;
					}
				}
				/*std::cout << "\nTIME COUNT: " << time_count << std::endl;
				for (int j=0; j<N_nodes_const; j++) {
					for (int i=0; i<N_nodes_const; i++) {
						std::cout << the_matrix[j][i] << "\t";
					}
					std::cout << std::endl;
				}*/

				//do backward substitution
				temperature_spatial[N_nodes_const - 1] = temperature_final;
				for (int i = N_nodes_const - 2; i>0; i--) {
					temperature_spatial[i] = (temperature_time[time_count][i] - the_matrix[i][i+1]*temperature_spatial[i+1])/the_matrix[i][i];
				}



				time_count = time_count + 1;
				std::cout << "step[" << time_count<< "]: "; 
				for (int i = 0; i<N_nodes_const; i++) {
					std::cout <<std::fixed<<std::setprecision(4) << temperature_spatial[i] << "\t";
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
				//std::cout << "ERRO: " <<error_computation << std::endl;
	
				if (time_count > 40) {
					std::cout << "GAGAL" << std::endl;
					break;
				}


			} while (error_computation > error_max);

		return temperature_spatial;


		}
};


int main() {

	Parameters par_1;
	//define input of par_1 object
	par_1.alpha			=	0.1;
	par_1.delta_t			=	0.05;
	par_1.delta_x			=	0.1;
	par_1.error_max			=	0.00001;
	par_1.temperature_initial	=	100;
	par_1.temperature_final		=	300;
	par_1.N_iter_max		=	100;
	par_1.N_nodes			=	20;

	TDMA_SOLVER_DIRICHLET tdma_1;

	std::vector<double> first = tdma_1.temperature_computation(par_1.alpha, par_1.delta_t, par_1.delta_x, par_1.error_max, par_1.temperature_initial, par_1.temperature_final, par_1.N_iter_max, par_1.N_nodes);

//	for (int i = 0; i<par_1.N_nodes; i++) {
//		std::cout << first[i] << "\t\t" ;
//	
//	}
	std::cout << std::endl;


}
