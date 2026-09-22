
  //Mito formulas
  //Uniporter - Song supplement 3.1 
  static constexpr double Bm_mito = 0.01 ; //Assume MPTP closed
  static constexpr double RTzF = (8.315 * 308) / (96.5 * 2); //z for Ca of +2  
  static constexpr double po_mito = 0.9 ; //MCU open probability 
  static constexpr int N_mcu = 200 ; //Number of MCUs per mito 
  static constexpr double gMCU_max = 8.1 ; //pS
  static constexpr double Km_uni = 19000 ; //micromolar
  static constexpr double zFvmyo = 2.0 * 96.5 * 18 ; //z of 2 for Ca, 18 microl for vmyo

  //Mito NCX - Song supplement 3.2 
  static constexpr double vNCX_max = 0.0035 ; //microM * ms^-1
  static constexpr double mito_v0 = 91 ; //mV
  static constexpr double kNa = 9.4 ; //mM  
  static constexpr double mito_ncx_n = 3 ; //NCX cooperativity for Nai
  static constexpr double kCa = 0.375 ; //micromolar 
  static constexpr double bfRT = (0.5 * 96.5) / (8.316 * 308); //constant b of 0.5 

  //Mito Membrane Potential (psi) - Song Supplement 3.6
  static constexpr double V_mitos = 3.5 ; //mV*ms*-1
  static constexpr double C_mito = 1.812 ; //microMolar * mv-1
  static constexpr double k_mitou = 0.0192 ; //ms-1

  //ATP constants - Song supplement 4.1 
  static constexpr double TAN = 7000 ; //Total nucleotide content
  static constexpr double ADP_buffer_rate = 0.025; //Proportion of free cyto ADP 
  static constexpr double k_ATPconsum = 0.01 ; 
  static constexpr double k_dpspconsum = 200 ; 
  static constexpr double k_ATPase = 0.16 ; 
  static constexpr double kd_deltam = 150;
  static constexpr double kd_psn = 0.02 ; 
  static constexpr double kd_psp = 0.5 ; 
  static constexpr double DATP = 0.25 ; //diffusion coefficient
  static constexpr double dx = 1.84 ; //distance of diffusion