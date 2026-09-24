// #ifdef _OPENMP
#include <omp.h>
// #endif

#include <iostream>
#include <fstream>
using namespace std;
#include "subcell.hpp"
#include "SAN_elecphysio.hpp"
#include "input_output.h"
// #include <direct.h>

// #define CAV3_INCX
// #define CAV3_ICaT
// #define CAV3_ICaL
// #define RAND_SEED
// #define OUT_CJ

int main(int argc, char *argv[]) {


	int LTCC_alpha = 8;
	int LTCC_beta = 0;
	int LTCC_gamma = 0;
	// srand(50);
	CSubcell sc;
	sc.finemesh = 1;
	sc.nx = 34;
	sc.ny = 8;
	sc.nz = 8;
	double dt = 0.01;
	sc.setdt(dt);


	bool Na_clamp = true;



	SAN_elecphysio Cell;

	int Tn = 31000.0 / dt;
	ofstream os("ci.txt");
	double v = -80;

#ifdef CAV3_INCX
	if (argc == 3) {

		sc.NCXalpha = atof(argv[1]);   // propotion of NCX in the cleft area
		double NCX_beta = atof(argv[2]);  // propotion of NCX in the subsarcolemal area

		sc.NCX_gamma = 1 - (sc.NCXalpha + NCX_beta);

		if (sc.NCX_gamma  > 1 or sc.NCX_gamma  < 0) {
			std::cerr << " incorrect setting of NCX alpha and beta parameters!!!!" << std::endl;
			std::cerr << " incorrect setting of NCX alpha and beta parameters!!!!" << std::endl;
			std::cerr << " incorrect setting of NCX alpha and beta parameters!!!!" << std::endl;
			std::exit(0);
		}
	}
#endif



#ifdef CAV3_ICaT
	if (argc == 3) {

		sc.alpha_CaT = atof(argv[1]);   // propotion of ICaT in the cleft area
		double CaT_beta = atof(argv[2]);  // propotion of ICaT in the subsarcolemal area

		sc.gamma_CaT = 1 - (sc.alpha_CaT + CaT_beta);

		if (sc.gamma_CaT  > 1 or sc.gamma_CaT  < 0) {
			std::cerr << " incorrect setting of ICAT alpha and beta parameters!!!!" << std::endl;
			std::cerr << " incorrect setting of ICAT alpha and beta parameters!!!!" << std::endl;
			std::cerr << " incorrect setting of ICAT alpha and beta parameters!!!!" << std::endl;
			std::exit(0);
		}
	}
#endif



#ifdef CAV3_ICaL


	if (argc == 3) {

		LTCC_alpha = atoi(argv[1]);   // propotion of ICaT in the cleft area
		LTCC_beta = atoi(argv[2]);  // propotion of ICaT in the subsarcolemal area

		LTCC_gamma = 8 - (LTCC_alpha + LTCC_beta); 

		if (LTCC_gamma  > 8 or LTCC_gamma  < 0) {
			std::cerr << " incorrect setting of LTCC alpha and beta parameters!!!!" << std::endl;
			std::cerr << " incorrect setting of LTCC alpha and beta parameters!!!!" << std::endl;
			std::cerr << " incorrect setting of LTCC alpha and beta parameters!!!!" << std::endl;
			std::exit(0);
		}
	}
#endif

	int random_seed = 0;

#ifdef RAND_SEED
	
	random_seed = atoi(argv[1]);

#endif



	sc.CLAMP_Cai = false;  
	sc.init(0.1, 800, LTCC_alpha, LTCC_gamma,random_seed);

	std::cout << sc.NUM_Ttubule << std::endl; 

	//Initialize CRU properties 
	sc.set_CRU_type();
	sc.assign_mito(); 
	sc.assign_producer_CRU() ;

	sc.output_map(sc.CRU_type, "CRU_type.vtk", 0, 0, 0);
	sc.output_map(sc.tubule_flag, "tubule_flag.vtk", 0, 0, 0);
	sc.output_map(sc.CRU_mito_assignment, "mito_idx.vtk", 0, 0, 0);
	sc.output_map(sc.CRU_producer_status, "prod_idx.vtk", sc.nx_mito, sc.ny_mito, sc.nz_mito);

	// //Temporary loop to print mito indeces 
	// for (int i = 0; i < sc.n; ++i){
	// 	std::cout << sc.CRU_mito_assignment[i] << "\n" ;
	// 	}
	//Testing averages across prod crus
	std::ofstream mito_summary("mito_summary.txt");
	mito_summary << "time\tavg_ATP\tavg_ADP\tavg_ca_mito\tavg_psi_mito\n";

	std::ofstream mito0_trace("mito0_trace.txt");
	std::ofstream atp_linescan("atp_linescan.txt");
	mito0_trace << "time\tcp0\tJ_uni0\tjNCX_m0\tca_mito0\tpsi_mito0\n";
	// print CRu type with producer mito
	int count_type0 = 0, count_type1 = 0, count_type2 = 0;
	for (int id_mito = 0; id_mito < sc.n_mito; ++id_mito) {
		int producer = sc.CRU_producer_status[id_mito];
		int type = sc.CRU_type[producer];
		if (type == 0) count_type0++;
		else if (type == 1) count_type1++;
		else if (type == 2) count_type2++;
	}

	std::cout << "type 0 (interior): " << count_type0 << std::endl;
	std::cout << "type 1 (near-boundary): " << count_type1 << std::endl;
	std::cout << "type 2 (boundary): " << count_type2 << std::endl;
	std::cout << "total: " << (count_type0 + count_type1 + count_type2) << std::endl;
	
	//Visual for CRUs attached to mito vs whole grid
	int *producer_highlight = new int[sc.n];
	for (int i = 0; i < sc.n; ++i) producer_highlight[i] = 0;

	for (int id_mito = 0; id_mito < sc.n_mito; ++id_mito) {
	producer_highlight[sc.CRU_producer_status[id_mito]] = 1;
	}

	sc.output_map(producer_highlight, "producer_highlight.vtk", 0,0,0);  // default CRU-grid dims
	delete [] producer_highlight;
	// to simulate ion current blockade
	// sc.ncx_scale = 0.4;
	// sc.ICaT_scale = 0.4;
	// sc.ICaL_scale = 0.4;
	// Cell.If_scale = 0.4; // If parameters in cell instead of sc.

	

	for (int tn = 0; tn < Tn; tn++)
	{

		// operator splitting pt 1


		double t = tn * dt;

		// note that Cm = 0.025 nF from SAN_elecphysio.hpp // 16:21:14, Mon, 04-May-2020, By Haibo
		Cell.update_Na_and_K_currents(t);


		// time capacitance of the cell here, capacitance = 0.025 in the original Kharche model; 
		Cell.icat = sc.ICaT_tot * 0.025;

		Cell.ical13 = sc.ica_stan * 0.025;
		Cell.ibca = sc.icabk_stan * 0.025;
		Cell.inaca = sc.incx_stan * 0.025;
		Cell.icap = sc.ipca_stan * 0.025;

		if (not Na_clamp) {
			Cell.update_Na_K_concentration(t);  // update Nai here as well. // 17:04:27, Wed, 29-April-2020, By Haibo
		}
		Cell.com_total_current(t);

		// Cell.update(t);
		Cell.update_state_FE(dt / 2.0);
		/*if (t > 22900 and t < 23900)
			Cell.y[36] = -65;*/
		v = Cell.y[36];

		// i_CaT is implemented in the pace function, 
		// solved with a single dt here.
		if (Na_clamp) {

			sc.pace(v, 10.0);  // fix Na concentration
		} else {
			sc.pace(v, Cell.y[34]);
		}


		// operator splitting pt 2

		Cell.icat = sc.ICaT_tot * 0.025;

		// Cell.update(t+dt/2.0);
		Cell.update_Na_and_K_currents(t);

		// Cell.update_Ca_currents(t);
		// Cell.com_INaCa(t);
		Cell.ical13 = sc.ica_stan * 0.025;
		Cell.ibca = sc.icabk_stan * 0.025;
		Cell.inaca = sc.incx_stan * 0.025;
		Cell.icap = sc.ipca_stan * 0.025;

		if (not Na_clamp) {
			Cell.update_Na_K_concentration(t);  // update Nai here as well. // 17:04:27, Wed, 29-April-2020, By Haibo
		}
		Cell.com_total_current(t);
		Cell.update_state_FE(dt / 2.0);

		/*if (t > 22900 and t < 23900)
			Cell.y[36] = -65;*/


		if (tn % 100 == 0)
		{
			cout << t << "\t" << sc.ci[0]  << "\t" << sc.cp[0] << endl;
			// os << t << "\t" << v << "\t" << sc.ica_stan << "\t" << sc.num_open_ICaL << "\t" << Cell.ical12 / 0.025 << "\t" << Cell.ical13 / 0.025 << "\t" ;
			double avg_ADP = ADP_buffer_rate * (TAN - sc.avg_ATP);
			mito_summary << t << "\t" << sc.avg_ATP << "\t" << avg_ADP
						<< "\t" << sc.avg_ca_mito << "\t" << sc.avg_psi_mito << "\n";
			mito_summary.flush();
			mito0_trace << t << "\t" << sc.trace_cp0 << "\t" << sc.trace_Juni0 << "\t"
					<< sc.trace_jncx0 << "\t" << sc.trace_ca_mito0 << "\t" << sc.trace_psi_mito0 << "\n";
			mito0_trace.flush();
			int j_mid = sc.ny / 2;
			int k_mid = sc.nz / 2;
			int i_fixed = sc.nx / 2;   // x index, arbitrary (all x are producer-level)
			int k_fixed = sc.nz / 2;   // z index — must be even (producer z-plane)
			for (int j = 0; j < sc.ny; ++j) {
				atp_linescan << sc.ATP_cyto[i_fixed + j*sc.nx + k_fixed*(sc.nx*sc.ny)];
				if (j < sc.ny - 1) atp_linescan << "\t";
			}
			atp_linescan << "\n";
			atp_linescan.flush();
			os << t << " " << Cell.y[37 - 1] << " " << sc.compute_avg_ci()  << " " << sc.compute_avg_cnsr()
			   << " " << Cell.ih / 0.025   //5
			   << " " << Cell.ina_ttxs / 0.025
			   << " " << Cell.ina_ttxr / 0.025
			   << " " << Cell.ical12 / 0.025
			   << " " << Cell.ical13 / 0.025
			   << " " << Cell.iks / 0.025 //10
			   << " " << Cell.ikr / 0.025
			   << " " << Cell.ik1 / 0.025
			   << " " << Cell.ist / 0.025  // 13
			   << " " << Cell.ib / 0.025
			   << " " << Cell.icat / 0.025
			   << " " << Cell.inak / 0.025 //16
			   << " " << Cell.isus / 0.025
			   << " " << Cell.inaca / 0.025
			   << " " << Cell.ito / 0.025
			   << " " <<  Cell.ibna / 0.025  // 20
			   << " " <<  Cell.ibca / 0.025 // 21
			   << " " <<  Cell.ibk / 0.025  // 22
			   << " " << Cell.icap / 0.025
			   << " " << sc.ica_stan   // 24
			   << " " << sc.incx_stan
			   << " " << sc.icabk_stan
			   << " " << sc.ipca_stan
			   << " " << sc.compute_avg_cp()
			   << " " << sc.compute_avg_cs()  // 29
			   << " " << Cell.y[34]
			   << " " <<  sc.compute_avg_cjsr()
			   // << " " << sc.incx_stan
			   << std::endl;
			// os << sc.compute_avg_ci() << "\t" << sc.compute_avg_cnsr() << endl;



			/* for (int id=0;id<sc.nx;id++) {
			   os <<sc.ci[id]<<"\t";
			 }*/
			/*     for (int id_nx = 0; id_nx < sc.nx; id_nx++)
			     {
			       int id = id_nx + 5 * sc.nx + 5 * sc.nx * sc.ny;
			       os << sc.ci[id] << "\t";
			     }
			     for (int in_ny = 0; in_ny < sc.nx; in_ny++)
			     {
			       int id = 5 + in_ny * sc.nx + 5 * sc.nx * sc.ny;
			       os << sc.ci[id] << "\t";
			     }

			     for (int id_nx = 0; id_nx < sc.nx; id_nx++)
			     {
			       int id = id_nx + 5 * sc.nx + 5 * sc.nx * sc.ny;
			       os << sc.cnsr[id] << "\t";
			     }
			*/
			if (t > 1000 and tn % 100 == 0) {

				char filename[1000];
				sprintf(filename, "Data/ci_%04d.bin", tn / 100);

				ofstream out_ci(filename);

				#ifdef OUT_CJ
				sprintf(filename, "Data/cj_%04d.bin", tn / 100);
				ofstream out_cj(filename);
				#endif
				if ( not out_ci.is_open()) {
					std::cerr << filename << " NOT opened !!! making a new one now!!!" << std::endl;
					// std::exit(0);
					// mkdir("Data/");
					const int dir_err = system("mkdir -p Data");
				}
				for (int k = 0; k < sc.nz; ++k)
					for (int i = 0; i < sc.ny; ++i)
					{
						for (int j = 0; j < sc.nx; ++j)
						{
							int id = j + i * sc.nx + k * sc.nx * sc.ny;
							#ifdef OUT_CJ
							out_cj << sc.cp[id] << "\t";
							#endif

							out_ci << sc.ci[id] << "\t";
						}
							#ifdef OUT_CJ
						out_cj << std::endl;
							#endif

						out_ci << std::endl;

					}
				#ifdef OUT_CJ
				out_cj.close();
				#endif

				out_ci.close();

				// output_double_array_bin(filename, sc.ci, sc.nx * sc.nx * sc.ny);
			}
		}
	}
	sc.output_map(sc.ATP_cyto, "atp_cyto_final.vtk", 0, 0, 0); 
	sc.output_map(sc.ATP_prod_rate, "atp_prod_rate_final.vtk", 0, 0, 0);
	mito_summary.close() ;
	atp_linescan.close();
	return 0;
}
