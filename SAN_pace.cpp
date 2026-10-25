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

int main(int argc, char *argv[])
{

	int LTCC_alpha = 8;
	int LTCC_beta = 0;
	int LTCC_gamma = 0;
	// srand(50);
	CSubcell sc;
	sc.finemesh = 1;
	sc.nx = 34; //Original
	//sc.nx = 36; //Upsize to test different CRU:Mito ratios (Necessary if testing 6,8,12,16, or 24:1)
	sc.ny = 8;  //Original
	sc.nz = 8;  //Original
	double dt = 0.01;
	sc.setdt(dt);

	//Adjust boolean clamps
	bool Na_clamp = true;
	//AP clamp with pregenerated waveforms 
	bool AP_clamp_on = false;  
	double AP_clamp_CL = 5000.0 ;  //106 (validation), 135, 150, 190, 235, 250, 285, 350 ;
	//Must  extend simulation time based on CL for AP clamp to match number of beats for a 100 beat simulation with BCL 

	//Mito Grid Dimensions - Ratios - 4:1 (1x2x2) ; 6:1 (3x2x1) ; 8:1 (2x2x2) ; 12:1 (3x2x2) : 16:1 (1x4x4) ; 24:1 (6x2x2)
	sc.nxover_nxmito = 1; 
	sc.nyover_nymito = 2 ;
	sc.nzover_nzmito = 2 ;

	//Clamp ATP for whole simulation
	sc.ATP_clamp_on = false; 
	sc.ATP_clamp = 20; 
	//Consistently impair ATP production by scaling factor for every producer CRU 
	sc.ATP_impair_on = false ;
	sc.ATP_impair_scale = 0.2 ; 

	SAN_elecphysio Cell;
	int Tn = 60000.0 / dt;
	ofstream os("ci.txt");
	double v = -80; 
	
#ifdef CAV3_INCX
	if (argc == 3)
	{

		sc.NCXalpha = atof(argv[1]);	 // propotion of NCX in the cleft area
		double NCX_beta = atof(argv[2]); // propotion of NCX in the subsarcolemal area

		sc.NCX_gamma = 1 - (sc.NCXalpha + NCX_beta);

		if (sc.NCX_gamma > 1 or sc.NCX_gamma < 0)
		{
			std::cerr << " incorrect setting of NCX alpha and beta parameters!!!!" << std::endl;
			std::cerr << " incorrect setting of NCX alpha and beta parameters!!!!" << std::endl;
			std::cerr << " incorrect setting of NCX alpha and beta parameters!!!!" << std::endl;
			std::exit(0);
		}
	}
#endif

#ifdef CAV3_ICaT
	if (argc == 3)
	{

		sc.alpha_CaT = atof(argv[1]);	 // propotion of ICaT in the cleft area
		double CaT_beta = atof(argv[2]); // propotion of ICaT in the subsarcolemal area

		sc.gamma_CaT = 1 - (sc.alpha_CaT + CaT_beta);

		if (sc.gamma_CaT > 1 or sc.gamma_CaT < 0)
		{
			std::cerr << " incorrect setting of ICAT alpha and beta parameters!!!!" << std::endl;
			std::cerr << " incorrect setting of ICAT alpha and beta parameters!!!!" << std::endl;
			std::cerr << " incorrect setting of ICAT alpha and beta parameters!!!!" << std::endl;
			std::exit(0);
		}
	}
#endif

#ifdef CAV3_ICaL

	if (argc == 3)
	{

		LTCC_alpha = atoi(argv[1]); // propotion of ICaT in the cleft area
		LTCC_beta = atoi(argv[2]);	// propotion of ICaT in the subsarcolemal area

		LTCC_gamma = 8 - (LTCC_alpha + LTCC_beta);

		if (LTCC_gamma > 8 or LTCC_gamma < 0)
		{
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
	//sc.setlayer(4);   // TEST: flag all CRUs as peripheral
	//Can override mito constructor defaults here if desired, must be before sc.init
	sc.mito = Mito(4,4,1,1) ; //MCU_scale, vNCX scale, MCU_block_activity, mNCX_block_activity defaults are 4,4,1,1 - Anthony 10-25-2026
	//For example, simulating 50% MCU block (4,4,0.5,1), 50% mNCX block: (4,4,1,0.5) ;
	sc.init(0.1, 800, LTCC_alpha, LTCC_gamma, random_seed);

	std::cout << sc.NUM_Ttubule << std::endl;
	std::cout << "  vNCX_scale " << sc.mito.vNCX_scale << "kMCU scale " << sc.mito.Km_uni_scale << "MCU_block_scale " 
	<< sc.mito.MCU_block_scale << "mNCX block scale " << sc.mito.mNCX_block_scale << std::endl;
	// Initialize CRU properties
	sc.set_CRU_type();
	sc.assign_mito();
	sc.assign_producer_CRU();

	sc.output_map(sc.CRU_type, "CRU_type.vtk", 0, 0, 0);
	sc.output_map(sc.tubule_flag, "tubule_flag.vtk", 0, 0, 0);
	sc.output_map(sc.CRU_mito_assignment, "mito_idx.vtk", 0, 0, 0);
	sc.output_map(sc.CRU_producer_status, "prod_idx.vtk", sc.nx_mito, sc.ny_mito, sc.nz_mito);
	sc.output_map(sc.ATP_impair, "atp_impair_map.vtk", sc.nx_mito, sc.ny_mito, sc.nz_mito) ;
	sc.output_map(sc.impair_status_map, "producer_status_map.vtk", 0, 0, 0);

	// //Temporary loop to print mito indeces
	// for (int i = 0; i < sc.n; ++i){
	// 	std::cout << sc.CRU_mito_assignment[i] << "\n" ;
	// 	}
	// Testing averages across prod crus
	std::ofstream mito_summary("mito_summary.txt");
	std::ofstream mito0_trace("mito0_trace.txt");
	std::ofstream atp_group_trace("atp_group_trace.txt");
	std::ofstream mito_type("mito_by_type.txt");

	mito_summary << "time\tavg_ATP\tavg_ADP\tavg_ca_mito\tavg_psi_mito\tavg_cai\tavg_cp_prod\tavg_J_uni\tavg_jNCX_m\tavg_Jflux_cp\tavg_NCXflux_ci\n";
	mito0_trace << "time\tcp0\tcai0\tJ_uni0\tjNCX_m0\tca_mito0\tpsi_mito0\tatp0\tadp0\tprod0\tconsum0\tdiff0\tmitoJfluxCp0\tmitoNCXfluxCi0\n";
	atp_group_trace << "time\tATP_nonprod\tATP_prod_unimpaired\tATP_prod_impaired\n"; // header, once before loop
	mito_type << "time\tca_mito_0\tca_mito_1\tca_mito_2\tpsi_0\tpsi_1\tpsi_2\tcp_0\tcp_1\tcp_2\tJuni_0\tJuni_1\tJuni_2\n" ;
	std::ofstream atp_linescan("atp_linescan.txt");
	std::ofstream atp_full("atp_full_grid.txt");
	std::ofstream ikatp_inak_serca_trace("ikatp_inak_serca_trace.txt");
	//std::ofstream vm_hr("vm_hr.txt");
	ikatp_inak_serca_trace << "time\tATP_ave\tp_kATP\tikATP\tinak\tIup_avg\tnai\n";

	// print CRu type with producer mito
	int count_type0 = 0, count_type1 = 0, count_type2 = 0;
	for (int id_mito = 0; id_mito < sc.n_mito; ++id_mito)
	{
		int producer = sc.CRU_producer_status[id_mito];
		int type = sc.CRU_type[producer];
		if (type == 0)
			count_type0++;
		else if (type == 1)
			count_type1++;
		else if (type == 2)
			count_type2++;
	}

	std::cout << "type 0 (interior): " << count_type0 << std::endl;
	std::cout << "type 1 (near-boundary): " << count_type1 << std::endl;
	std::cout << "type 2 (boundary): " << count_type2 << std::endl;
	std::cout << "total: " << (count_type0 + count_type1 + count_type2) << std::endl;

	// Visual for CRUs attached to mito vs whole grid
	int *producer_highlight = new int[sc.n];
	for (int i = 0; i < sc.n; ++i)
		producer_highlight[i] = 0;

	for (int id_mito = 0; id_mito < sc.n_mito; ++id_mito)
	{
		producer_highlight[sc.CRU_producer_status[id_mito]] = 1;
	}

	sc.output_map(producer_highlight, "producer_highlight.vtk", 0, 0, 0); // default CRU-grid dims
	std::ofstream producer_mask("producer_mask.txt");
	for (int i = 0; i < sc.n; ++i)
	{
		producer_mask << producer_highlight[i];
		if (i < sc.n - 1)
			producer_mask << "\t";
	}
	producer_mask << "\n";
	producer_mask.close();

	// [added] Spatial snapshot output: per-mitochondrion Ca, Psi and producer-CRU cleft Ca,
	// written every 1 ms for the last snap_window_ms of the run (for linescans / heatmaps).
	// mito_map.txt gives, for each mito index, its grid position and producer CRU type.
	const double snap_window_ms = 2000.0;
	std::ofstream mito_map("mito_map.txt");
	mito_map << "id_mito\tmx\tmy\tmz\tprod_id\tCRU_type\n";
	for (int id_mito = 0; id_mito < sc.n_mito; ++id_mito)
	{
		int mx = id_mito % sc.nx_mito;
		int my = (id_mito / sc.nx_mito) % sc.ny_mito;
		int mz = id_mito / (sc.nx_mito * sc.ny_mito);
		int pid = sc.CRU_producer_status[id_mito];
		mito_map << id_mito << "\t" << mx << "\t" << my << "\t" << mz << "\t" << pid
				 << "\t" << sc.CRU_type[pid] << "\n";
	}
	mito_map.close();
	std::ofstream mito_ca_grid("mito_ca_grid.txt");     // row: t, ca_mito[0..n_mito-1]   (uM)
	std::ofstream mito_psi_grid("mito_psi_grid.txt");   // row: t, psi_mito[0..n_mito-1]  (mV)
	std::ofstream mito_cp_grid("mito_cp_grid.txt");     // row: t, cp[producer of each mito] (uM)

	delete[] producer_highlight;

	// to simulate ion current blockade
	//sc.ncx_scale = 0.85;
	//sc.ICaT_scale = 0.4;
	//sc.ICaL_scale = 0.4;
	//Cell.If_scale = 0.4; // If parameters in cell instead of sc.
	//Cell.ist_scale = 0.4 ; //ist pars in cell as well 

	//load appropriate AP clamp waveform if defined 
	std::vector<double> AP_set;
	if (AP_clamp_on) {
		std::string ap_file = "pool_AP/AP_" + std::to_string((int)AP_clamp_CL) + ".txt";
		std::ifstream inAP(ap_file);
		double tmp;
		while (inAP >> tmp) AP_set.push_back(tmp);		

		if (AP_set.empty()) {
			std::cerr << "AP file not found or empty: " << ap_file << std::endl;
			return 1;
		}
		if (AP_set.size() != (size_t)std::lround(AP_clamp_CL / dt)) {
			std::cerr << "AP file length " << AP_set.size()
			<< " != CL/dt " << std::lround(AP_clamp_CL / dt) << std::endl;
			return 1;
		}
	}

	//Integration loop
	for (int tn = 0; tn < Tn; tn++)
	{

		// operator splitting pt 1

		double t = tn * dt;
		//AP clamp conditions
		double v_clamp = 0.0 ;
		if(AP_clamp_on) {
			v_clamp = AP_set[tn % AP_set.size()] ;
			Cell.y[36] = v_clamp ;
		}
		// note that Cm = 0.025 nF from SAN_elecphysio.hpp // 16:21:14, Mon, 04-May-2020, By Haibo
		Cell.update_Na_and_K_currents(t, sc.avg_ATP, sc.avg_ADP);

		// time capacitance of the cell here, capacitance = 0.025 in the original Kharche model;
		Cell.icat = sc.ICaT_tot * 0.025;

		Cell.ical13 = sc.ica_stan * 0.025;
		Cell.ibca = sc.icabk_stan * 0.025;
		Cell.inaca = sc.incx_stan * 0.025;
		Cell.icap = sc.ipca_stan * 0.025;

		if (not Na_clamp)
		{
			Cell.update_Na_K_concentration(t); // update Nai here as well. // 17:04:27, Wed, 29-April-2020, By Haibo
		}
		Cell.com_total_current(t);

		// Cell.update(t);
		Cell.update_state_FE(dt / 2.0);
		/*if (t > 22900 and t < 23900)
			Cell.y[36] = -65;*/
		if(AP_clamp_on) Cell.y[36] = v_clamp; //undo Vm integration 
		v = Cell.y[36];
		//if (!AP_clamp_on && t >= Tn * dt - 2000.0 && tn % 10 == 0)
			//vm_hr << t << "\t" << v << "\n";
		// i_CaT is implemented in the pace function,
		// solved with a single dt here.
		if (Na_clamp)
		{

			sc.pace(v, 10.0); // fix Na concentration
		}
		else
		{
			sc.pace(v, Cell.y[34]);
		}

		// operator splitting pt 2

		Cell.icat = sc.ICaT_tot * 0.025;

		// Cell.update(t+dt/2.0);
		Cell.update_Na_and_K_currents(t, sc.avg_ATP, sc.avg_ADP);

		// Cell.update_Ca_currents(t);
		// Cell.com_INaCa(t);
		Cell.ical13 = sc.ica_stan * 0.025;
		Cell.ibca = sc.icabk_stan * 0.025;
		Cell.inaca = sc.incx_stan * 0.025;
		Cell.icap = sc.ipca_stan * 0.025;

		if (not Na_clamp)
		{
			Cell.update_Na_K_concentration(t); // update Nai here as well. // 17:04:27, Wed, 29-April-2020, By Haibo
		}
		Cell.com_total_current(t);
		Cell.update_state_FE(dt / 2.0);
		if (AP_clamp_on) Cell.y[36] = v_clamp;   // reset 
		v = Cell.y[36];

		/*if (t > 22900 and t < 23900)
			Cell.y[36] = -65;*/

		if (tn % 100 == 0)
		//if (tn % 1000 == 0) //for timestep = 0.001
		{
			// [added] spatial snapshots, last snap_window_ms only
			if (t >= Tn * dt - snap_window_ms)
			{
				mito_ca_grid << t;  mito_psi_grid << t;  mito_cp_grid << t;
				for (int id_mito = 0; id_mito < sc.n_mito; ++id_mito)
				{
					mito_ca_grid  << "\t" << sc.ca_mito[id_mito];
					mito_psi_grid << "\t" << sc.psi_mito[id_mito];
					mito_cp_grid  << "\t" << sc.cp[sc.CRU_producer_status[id_mito]];
				}
				mito_ca_grid << "\n";  mito_psi_grid << "\n";  mito_cp_grid << "\n";
			}
			cout << t << "\t" << sc.ci[0] << "\t" << sc.cp[0] << endl;
			// os << t << "\t" << v << "\t" << sc.ica_stan << "\t" << sc.num_open_ICaL << "\t" << Cell.ical12 / 0.025 << "\t" << Cell.ical13 / 0.025 << "\t" ;
			double avg_ADP = sc.mito.ADP_buffer_rate * (sc.mito.TAN - sc.avg_ATP);

			mito_summary << t << "\t" << sc.avg_ATP << "\t" << avg_ADP
						 << "\t" << sc.avg_ca_mito << "\t" << sc.avg_psi_mito
						 << "\t" << sc.compute_avg_ci() << "\t" << sc.avg_cp_prod
						 << "\t" << sc.avg_J_uni << "\t" << sc.avg_jNCX_m
						 << "\t" << sc.avg_Jflux_cp << "\t" << sc.avg_NCXflux_ci << "\n";
			mito_summary.flush();
			mito0_trace << t << "\t" << sc.trace_cp0 << "\t" << sc.trace_cai0 << "\t" << sc.trace_Juni0 << "\t"
							<< sc.trace_jncx0 << "\t" << sc.trace_ca_mito0 << "\t" << sc.trace_psi_mito0
							<< "\t" << sc.trace_atp0 << "\t" << sc.trace_adp0
							<< "\t" << sc.trace_prod0 << "\t" << sc.trace_consum0 << "\t" << sc.trace_diff0
							<< "\t" << sc.trace_mitoJfluxCp0 << "\t" << sc.trace_mitoNCXfluxCi0 << "\n";
			mito0_trace.flush();
			atp_group_trace << t << "\t" << sc.avg_ATP_nonprod << "\t" << sc.avg_ATP_prod_ok
							<< "\t" << sc.avg_ATP_prod_impaired << "\n";
			atp_group_trace.flush();
			mito_type << t;
				for (int k = 0; k < 3; ++k) mito_type << "\t" << sc.avg_ca_mito_type[k];
				for (int k = 0; k < 3; ++k) mito_type << "\t" << sc.avg_psi_mito_type[k];
				for (int k = 0; k < 3; ++k) mito_type << "\t" << sc.avg_cp_type[k];
				for (int k = 0; k < 3; ++k) mito_type << "\t" << sc.avg_Juni_type[k];
			mito_type << "\n";
			mito_type.flush();
			for (int id = 0; id < sc.n; ++id)
			{
				atp_full << sc.ATP_cyto[id];
				if (id < sc.n - 1)
					atp_full << "\t";
			}
			atp_full << "\n";
			atp_full.flush();
			ikatp_inak_serca_trace << t << "\t" << sc.avg_ATP << "\t" << Cell.p_kATP
								   << "\t" << Cell.ikATP << "\t" << Cell.inak
								   << "\t" << sc.iupave << "\t" << Cell.y[34] << "\n"; //Cell.y[34] is nai
			ikatp_inak_serca_trace.flush();
			int j_mid = sc.ny / 2;
			int k_mid = sc.nz / 2;
			int i_fixed = sc.nx / 2; // x index, arbitrary (all x are producer-level)
			int k_fixed = sc.nz / 2; // z index — must be even (producer z-plane)
			for (int j = 0; j < sc.ny; ++j)
			{
				atp_linescan << sc.ATP_cyto[i_fixed + j * sc.nx + k_fixed * (sc.nx * sc.ny)];
				if (j < sc.ny - 1)
					atp_linescan << "\t";
			}
			atp_linescan << "\n";
			atp_linescan.flush();
			os << t << " " << Cell.y[37 - 1] << " " << sc.compute_avg_ci() << " " << sc.compute_avg_cnsr()
			   << " " << Cell.ih / 0.025 // 5
			   << " " << Cell.ina_ttxs / 0.025
			   << " " << Cell.ina_ttxr / 0.025
			   << " " << Cell.ical12 / 0.025
			   << " " << Cell.ical13 / 0.025
			   << " " << Cell.iks / 0.025 // 10
			   << " " << Cell.ikr / 0.025
			   << " " << Cell.ik1 / 0.025
			   << " " << Cell.ist / 0.025 // 13
			   << " " << Cell.ib / 0.025
			   << " " << Cell.icat / 0.025
			   << " " << Cell.inak / 0.025 // 16
			   << " " << Cell.isus / 0.025
			   << " " << Cell.inaca / 0.025
			   << " " << Cell.ito / 0.025
			   << " " << Cell.ibna / 0.025 // 20
			   << " " << Cell.ibca / 0.025 // 21
			   << " " << Cell.ibk / 0.025  // 22
			   << " " << Cell.icap / 0.025
			   << " " << sc.ica_stan // 24
			   << " " << sc.incx_stan
			   << " " << sc.icabk_stan
			   << " " << sc.ipca_stan
			   << " " << sc.compute_avg_cp()
			   << " " << sc.compute_avg_cs() // 29
			   << " " << Cell.y[34]
			   << " " << sc.compute_avg_cjsr()
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
			if (t > 1000 and tn % 100 == 0)
			//if (t > 1000 and tn % 1000 == 0) //for timestep = 0.001
			{

				char filename[1000];
				sprintf(filename, "Data/ci_%04d.bin", tn / 100);
				//sprintf(filename, "Data/ci_%04d.bin", tn / 1000); //for timestep = 0.001

				ofstream out_ci(filename);

#ifdef OUT_CJ
				sprintf(filename, "Data/cj_%04d.bin", tn / 100);
				ofstream out_cj(filename);
#endif
				if (not out_ci.is_open())
				{
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
	mito_summary.close();
	atp_linescan.close();
	atp_full.close();
	mito_ca_grid.close();  mito_psi_grid.close();  mito_cp_grid.close();   // [added]
	return 0;
}