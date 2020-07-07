	dN_1 = ( - D * N_1 ) + N_1  * ( mu_max_1 * S_glu / ( K_mu_glu + S_glu ) )

	dN_2 = ( - D * N_2 ) + N_2  * ( mu_max_2 * S_glu / ( K_mu_glu + S_glu ) ) - (  omega_max * (C_extra * B_1)**n_omega / ( k_omega_B_1**n_omega + (C_extra *  B_1)**n_omega )  )  * N_2 

	dS_glu = ( D * ( S0_glu - S_glu ) ) - ( mu_max_1 * S_glu / ( K_mu_glu + S_glu ) ) * N_1 * C_OD / g_1  - ( mu_max_2 * S_glu / ( K_mu_glu + S_glu ) ) * N_2 * C_OD / g_2 

	dB_1 = ( - D * B_1 ) +  kB_max_1  * ( (C_extra * A_1)**n_A_B_1 / ( K_A_B_1**n_A_B_1 + (C_extra * A_1)**n_A_B_1 ) ) * N_1 * C_OD / C_extra  +  kB_max_1  * ( K_A_B_1**n_A_B_1 / ( K_A_B_1**n_A_B_1 + (C_extra * A_1)**n_A_B_1 ) ) * N_2 * C_OD / C_extra 

	dA_1 = ( - D * A_1 ) + kA_1 * N_1 * C_OD / C_extra

