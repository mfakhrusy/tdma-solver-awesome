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
			std::vector<double> temperature_temp(N_nodes_const);

			double B = B_constant_computation(alpha, delta_t, delta_x);
			double A = A_constant_computation(alpha, delta_t, delta_x);

			//solve the A x = y
			//A is the matrix, x is temperature_spatial, y is temperature final + s*BC on both edges.
			//initialize the vector by initial conditions 
			for (int i=0; i<N_nodes_const; i++) {
				temperature_spatial[i] = 0.;
			}

			//initialize the temperature_time
			for (int i=0; i<N_nodes_const; i++) {
				temperature_time[0][i] = temperature_spatial[i];
			}
			std::cout << std::endl;

			//build the y column
			temperature_temp[0] = temperature_initial;
			for  (int i=1; i<N_nodes_const-1; i++) {
				if (i == 1) {
					temperature_temp[i] = temperature_time[0][i] + A*temperature_initial;
				}
				else if (i == N_nodes_const-2) {
					temperature_temp[i] = temperature_time[0][i] + A*temperature_final;
				}
				else {
					temperature_temp[i] = temperature_time[0][i];
				}
			}
			temperature_temp[N_nodes_const-1] = temperature_final;

			std::cout << "step[0]: ";
			for (int i = 0; i<N_nodes_const; i++) std::cout << temperature_temp[i] << " "; std::cout << "\n";

			//make the initial matrix
			std::vector<std::vector<double>> the_matrix(N_nodes_const,std::vector<double>(N_nodes_const));

			for (int j=1; j<N_nodes_const-1; j++) {	//vertical -> rows
				for (int i=1; i<N_nodes_const-1; i++) {	//horizontal -> columns

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


			int time_count = 0;
			double error_computation = 0.;

			//looping process
			do {

				// ------------------------ do the gauss method ------------------
				for (int j = 2; j<N_nodes_const-1; j++) {
					double temp_const = the_matrix[j][j-1]/the_matrix[j-1][j-1];
					
					// first, make a bottom triangle matrix on A matrix
					for (int i=0; i<N_nodes_const; i++) {
						the_matrix[j][i] = the_matrix[j][i] - the_matrix[j-1][i]*temp_const;
					}

					//do the same operation on x column and y column
					temperature_temp[j] = temperature_temp[j] - temperature_temp[j-1]*temp_const;

				}

				//do backward substitution
				temperature_spatial[N_nodes_const - 1] 	= temperature_final;
				temperature_spatial[0]			= temperature_initial;
				temperature_spatial[N_nodes_const - 2]	= temperature_temp[N_nodes_const - 2]/the_matrix[N_nodes_const-2][N_nodes_const-2];
				for (int i = N_nodes_const - 3; i>0; i--) {
					//temperature_spatial[i] = (temperature_time[time_count][i] - the_matrix[i][i+1]*temperature_spatial[i+1])/the_matrix[i][i];
					temperature_spatial[i] = (temperature_temp[i] - the_matrix[i][i+1]*temperature_spatial[i+1])/the_matrix[i][i];
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


				// --------- BUILD NEW MATRIX -----------
				for (int j=1; j<N_nodes_const-1; j++) {	//vertical -> rows
					for (int i=1; i<N_nodes_const-1; i++) {	//horizontal -> columns
	
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

				//build the y column
				temperature_temp[0] = temperature_initial;
				for (int i=1; i<N_nodes_const-1; i++) {
		//			temperature_temp[i] = temperature_time[time_count-1][i];

					if (i == 1) {
						temperature_temp[i] = temperature_time[time_count][i] + A*temperature_initial;
					}
					else if (i == N_nodes_const-2) {
						temperature_temp[i] = temperature_time[time_count][i] + A*temperature_final;
					}
					else {
						temperature_temp[i] = temperature_time[time_count][i];
					}



				}
				temperature_temp[N_nodes_const-1] = temperature_final;
				
				//compute error
				for (int i = 0; i <N_nodes_const; i++) {
					double temp = abs(temperature_time[time_count][i] - temperature_time[time_count - 1][i]);
					error_computation = (error_computation + temp);		
				}
				error_computation = error_computation/N_nodes_const;
				//std::cout << "ERRO: " <<error_computation << std::endl;
	
				if (time_count > 200) {
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
	par_1.delta_x			=	0.05;
	par_1.error_max			=	0.00001;
	par_1.temperature_initial	=	10;
	par_1.temperature_final		=	300;
	par_1.N_iter_max		=	100;
	par_1.N_nodes			=	10;

	TDMA_SOLVER_DIRICHLET tdma_1;

	std::vector<double> first = tdma_1.temperature_computation(par_1.alpha, par_1.delta_t, par_1.delta_x, par_1.error_max, par_1.temperature_initial, par_1.temperature_final, par_1.N_iter_max, par_1.N_nodes);

//	for (int i = 0; i<par_1.N_nodes; i++) {
//		std::cout << first[i] << "\t\t" ;
//	
//	}
//	std::cout << std::endl;


}
