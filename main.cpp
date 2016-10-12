#include "global.hpp"

class TDMA_SOLVER {

	private:
		std::vector<double> temperature;	

	public:


};


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






	

}
