set	(SOURCE
	 test_integration.c
	 test_integration_convergence_support.c
	 test_integration_support_face.c
	)

set	(LIBS_DEPEND
	 Simulation
	 Test_Support_Containers
	)

set (LIB Test_Integration)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})


set (LIBS_BASE Test_Base Test_Support)

set (EXEC test_integration_mesh)
set (LIBS_DEPEND ${LIBS_BASE} General Mesh Simulation Test_Support_Containers)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "blended_2d_mixed.msh" "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_DG_ParametricMixed2D__ml0__p3")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "straight_2d_quad_periodic.msh" "euler/periodic_vortex/TEST_Euler_PeriodicVortex_QUAD__ml0__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "parametric_2d_quad_periodic_reflected.msh" "extern_mesh/TEST_parametric_2d_quad_periodic_reflected")

set (EXEC test_integration_fe_init)
set (LIBS_DEPEND ${LIBS_BASE} ${PETSC_LIBRARIES} Simulation Test_Support_Simulation)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "extern_mesh/TEST_blended_2d_mixed")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "extern_mesh/TEST_straight_2d_quad_periodic")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "extern_mesh/TEST_parametric_2d_quad_periodic_reflected")

set (EXEC test_integration_fe_init_refined)
set (LIBS_DEPEND ${LIBS_BASE} ${PETSC_LIBRARIES} Simulation Test_Integration Test_Support_Simulation)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "extern_mesh/TEST_h_refined_2d_tri_h1")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "extern_mesh/TEST_h_refined_2d_tri_h2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "extern_mesh/TEST_h_refined_2d_tri_h1+uniform")


set (EXEC test_integration_geometry)
set (LIBS_DEPEND ${LIBS_BASE} ${PETSC_LIBRARIES} Simulation Test_Support_Containers Test_Support_Simulation)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "extern_mesh/TEST_straight_2d_quad_periodic")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "extern_mesh/TEST_blended_2d_mixed")

set (EXEC test_integration_fluxes)
set (LIBS_DEPEND Test_Base Test_Integration)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "flux_advection_1d" "integration/TEST_Advection_Default_1d__ml0__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "flux_advection_2d" "integration/TEST_Advection_Default_2d__ml0__p2")
#add_test_DPG_w_path(${BIN_PATH_3D} ${EXEC} "flux_advection_3d" "integration/TEST_Advection_Default_3d__ml0__p2")
add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "flux_diffusion_1d" "integration/TEST_Diffusion_Default_1d__ml0__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "flux_diffusion_2d" "integration/TEST_Diffusion_Default_2d__ml0__p2")
#add_test_DPG_w_path(${BIN_PATH_3D} ${EXEC} "flux_diffusion_3d" "integration/TEST_Diffusion_Default_3d__ml0__p2")
#add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "flux_euler_1d" "integration/TEST_Euler_Default_1d__ml0__p0")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "flux_euler_2d" "integration/TEST_Euler_Default_2d__ml0__p0")
#add_test_DPG_w_path(${BIN_PATH_3D} ${EXEC} "flux_euler_3d" "integration/TEST_Euler_Default_3d__ml0__p0")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "flux_navier_stokes_2d" "integration/TEST_NavierStokes_Default_2d__ml0__p0") # TODO: Change this to use DPG such that Hessians are also tested.
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "flux_navier_stokes_sutherland_2d" "integration/TEST_NavierStokes_Sutherland_2d__ml0__p0") # TODO: Change this to use DPG such that Hessians are also tested.

set (EXEC test_integration_numerical_fluxes)
set (LIBS_DEPEND Test_Base Test_Integration Test_Support_Simulation)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "numerical_flux_advection_upwind_1d" "integration/numerical_fluxes/TEST_Advection_upwind_1d__ml1__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "numerical_flux_advection_upwind_2d" "integration/numerical_fluxes/TEST_Advection_upwind_2d__ml1__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "numerical_flux_euler_roe_pike_2d" "integration/numerical_fluxes/TEST_Euler_Roe_Pike_2d__ml1__p3")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "numerical_flux_euler_lax_friedrichs_2d" "integration/numerical_fluxes/TEST_Euler_Lax_Friedrichs_2d__ml1__p3")

set (EXEC test_integration_boundary_conditions)
set (LIBS_DEPEND Test_Base Test_Integration)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "boundary_conditions_advection_1d" "integration/TEST_Advection_Default_1d__ml0__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "boundary_conditions_advection_2d" "integration/TEST_Advection_Default_2d__ml0__p2")
#add_test_DPG_w_path(${BIN_PATH_3D} ${EXEC} "boundary_conditions_advection_3d" "integration/TEST_Advection_Default_3d__ml0__p2")
add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "boundary_conditions_diffusion_1d" "integration/TEST_Diffusion_Default_1d__ml0__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "boundary_conditions_diffusion_2d" "integration/TEST_Diffusion_Default_2d__ml0__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "boundary_conditions_euler_2d_riemann_slipwall" "integration/boundary_conditions/TEST_Euler_RiSl_2d__ml0__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "boundary_conditions_euler_2d_supersonic_in/out" "integration/boundary_conditions/TEST_Euler_SiSo_2d__ml0__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "boundary_conditions_euler_2d_back_p_total_tp" "integration/boundary_conditions/TEST_Euler_BpTtp_2d__ml0__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "boundary_conditions_navier_stokes_2d_noslip_general/adiabatic" "integration/boundary_conditions/TEST_NavierStokes_NSGA_2d__ml0__p2")

set (EXEC test_integration_non_conforming)
set (LIBS_DEPEND Test_Base Test_Integration)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "integration/TEST_non_conforming_mixed_2d")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "integration/TEST_non_conforming_unaligned_mixed_2d")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/free_stream/TEST_advection_free_stream_mixed2d_non_conforming_p_coarse")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "RESULTS_euler_joukowski_full_dg_2D_ar4_iso")

set (EXEC test_integration_free_stream)
set (LIBS_DEPEND ${LIBS_BASE} Core Simulation Test_Integration)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/free_stream/TEST_advection_free_stream_mixed2d__ml2__p3" "petsc_options_gmres_tol_1e-15")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/free_stream/TEST_advection_free_stream_mixed2d_unaligned__ml2__p3" "petsc_options_gmres_tol_1e-15")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/free_stream/TEST_advection_free_stream_mixed2d_trig_xl__ml2__p3" "petsc_options_gmres_tol_1e-15")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/free_stream/TEST_advection_free_stream_mixed2d_non_conforming" "petsc_options_gmres_tol_1e-15")
#add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/free_stream/TEST_advection_free_stream_mixed2d_non_conforming_p4_super" "petsc_options_gmres_tol_1e-15") # Expected to be failing, but with free-stream error converging optimally.

set (EXEC test_integration_restart_output)
set (LIBS_DEPEND ${LIBS_BASE} Core Simulation Test_Integration)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "integration/restart/TEST_Euler_SupersonicVortex_DG_restart_2D")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "integration/restart/TEST_euler_gaussian_bump_dg_restart_2d" "petsc_options_gmres_default")

set (EXEC test_integration_restart)
set (LIBS_DEPEND ${LIBS_BASE} Core Simulation Test_Integration)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "integration/restart/TEST_Euler_SupersonicVortex_DG_ParametricQuad2D")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "integration/restart/TEST_Euler_SupersonicVortex_DG_ParametricTri2D_p_le_2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "integration/restart/TEST_Euler_SupersonicVortex_DG_ParametricTri2D") # Failing (May 4th, 2018)
#add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "integration/restart/TEST_Euler_SupersonicVortex_DG_ParametricMixed2D") # Fix general TRI above then re-enable

# Potentially add test for equivalence between explicit/implicit computation of grad_coef.

set (EXEC test_integration_conservation)
set (LIBS_DEPEND ${LIBS_BASE} Core Simulation Test_Integration)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/default/dg/TEST_Advection_Default_DG_Mixed2D__ml2__p2" "petsc_options_gmres_default")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/default/dpg/TEST_Advection_Default_Conservative_DPG_Mixed2D__ml2__p2" "petsc_options_gmres_default")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_DG_ParametricMixed2D__ml0__p3" "petsc_options_gmres_default")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_Conservative_DPG_BlendedMixed2D__ml1__p2" "petsc_options_gmres_default")

set (EXEC test_integration_inf_sup)
set (LIBS_DEPEND ${LIBS_BASE} ${SLEPC_LIBRARIES} Simulation Test_Integration Test_Support_Simulation)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})


set (EXEC test_integration_linearization)
set (LIBS_DEPEND ${LIBS_BASE}
     Simulation
     Test_Integration
     Test_Support_General
     Test_Support_Simulation
    )
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "advection/default/dpg/TEST_Advection_Demkowicz_DPGII_LINE__ml2__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/peterson/dg/TEST_Advection_Peterson_TRI__ml2__p3")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/peterson/dpg/TEST_Advection_Peterson_DPG_TRI__ml1__p3")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/default/dpg/TEST_Advection_Default_Conservative_DPG_Mixed2D__ml2__p2")
add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "advection/default/opg/TEST_opg_advection_default__1d__ml4__p4")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/default/opg/TEST_opg_advection_default__mixed2d__ml2__p2")
add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "diffusion/steady/default/dg/TEST_Diffusion_Steady_Default_DG_LINE__ml2__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "diffusion/steady/default/dg/TEST_Diffusion_Steady_Default_DG_TRI__ml2__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "diffusion/steady/default/dg/TEST_Diffusion_Steady_Default_DG_QUAD__ml2__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "diffusion/steady/default/dg/TEST_Diffusion_Steady_Default_DG_Mixed2D__ml2__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_DG_ParametricMixed2D__ml0__p3")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_DPG_ParametricMixed2D__ml0__p3")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_DG_BlendedMixed2D__ml0__p3")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_Conservative_DPG_BlendedMixed2D__ml1__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/gaussian_bump/TEST_Euler_GaussianBump_ParametricMixed2D__ml0__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "navier_stokes/steady/taylor_couette/dg/TEST_NavierStokes_TaylorCouette_DG_ParametricMixed2D__ml0__p2")

set (EXEC test_integration_convergence)
set (LIBS_DEPEND ${LIBS_BASE} Core Simulation Test_Integration)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/peterson/dg/TEST_Advection_Peterson_TRI__ml0__p0" "petsc_options_gmres_default")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/peterson/dpg/TEST_Advection_Peterson_DPG_TRI__ml0__p0" "petsc_options_cg_ilu1") # Incomplete Cholesky was slower than LU.
add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "advection/default/dg/TEST_Advection_Default_DG_LINE" "petsc_options_gmres_default")
add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "advection/default/dpg/TEST_Advection_Demkowicz_DPGII_LINE" "petsc_options_cg_default")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/default/dpg/TEST_Advection_Default_DPG_TRI" "petsc_options_cg_ilu1")
add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "advection/default/opg/TEST_opg_advection_default__1d" "petsc_options_gmres_tol_1e-15")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/default/opgc/TEST_opgc_advection_default_xl__mixed2d" "petsc_options_gmres_tol_1e-15")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/free_stream/TEST_advection_free_stream_mixed2d_non_conforming_p4_super" "petsc_options_gmres_tol_1e-15")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/default/dpg/TEST_Advection_Default_DPG_QUAD" "petsc_options_cg_icc")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/default/dpg/TEST_Advection_Default_DPG_Mixed2D" "petsc_options_cg_icc")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/default/dpg/TEST_Advection_Default_Conservative_DPG_Mixed2D__ml0" "petsc_options_cg_default")
add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "diffusion/steady/default/dg/TEST_Diffusion_Steady_Default_DG_LINE" "petsc_options_cg_ilu1")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "diffusion/steady/default/dg/TEST_Diffusion_Steady_Default_DG_Mixed2D__ml0" "petsc_options_cg_ilu1")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/periodic_vortex/TEST_Euler_PeriodicVortex_QUAD__ml0__p2" "petsc_options_empty")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_DG_ParametricMixed2D" "petsc_options_gmres_default")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/NURBS_bump/TEST_Euler_NURBSbump_ParametricQUAD2D" "petsc_options_gmres_default")

#add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/vortex/TEST_opg_advection_vortex_mixed2d" "petsc_options_gmres_tol_1e-15")
#add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/NURBS/TEST_Euler_NURBS_ParametricQUAD2D" "petsc_options_gmres_default")
#add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/NURBS/TEST_Euler_NURBS_ParametricQUAD2D__ml1__p2" "petsc_options_gmres_default")
#add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_Conservative_DPG_BlendedMixed2D__ml1" "petsc_options_gmres_r120") # Having difficulty converging
#add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_DPG_ParametricMixed2D__ml0" "petsc_options_gmres_default")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_DG_BlendedMixed2D__ml1" "petsc_options_gmres_default")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/joukowski/TEST_Euler_Joukowski_DG_ParametricMixed2D" "petsc_options_gmres_default")
#add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/joukowski/TEST_Euler_Joukowski_Full_DG_ParametricMixed2D__ml1" "petsc_options_gmres_r120")
#add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_DPG_BlendedMixed2D__ml1" "petsc_options_gmres_default")
#add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/free_stream/TEST_Euler_FreeStream_DG_StraightTRI2D__ml0" "petsc_options_gmres_default")
#add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/free_stream/TEST_Euler_FreeStream_DG_ParametricMixed2D__ml0" "petsc_options_gmres_default")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "navier_stokes/steady/taylor_couette/dg/TEST_NavierStokes_TaylorCouette_DG_ParametricTRI" "petsc_options_gmres_default")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "navier_stokes/steady/taylor_couette/dg/TEST_NavierStokes_TaylorCouette_DG_ParametricQUAD" "petsc_options_gmres_default")

set (EXEC test_integration_adjoint)
set (LIBS_DEPEND ${LIBS_BASE} Core Simulation Test_Integration)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/NURBS_Airfoil/TEST_Euler_NURBSAirfoil_AdjointSensitivitiesTestsReference_ParametricQUAD2D__ml1__p2" "petsc_options_gmres_default")

set (EXEC test_integration_sensitivities)
set (LIBS_DEPEND ${LIBS_BASE} Core Simulation Test_Integration)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/NURBS_Airfoil/TEST_Euler_NURBSAirfoil_AdjointSensitivitiesTestsReference_ParametricQUAD2D__ml1__p2" "petsc_options_gmres_default")

set (EXEC test_integration_run)
set (LIBS_DEPEND ${LIBS_BASE} Core Simulation Test_Integration)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/NURBS_bump/TEST_Euler_NURBSbump_ParametricQUAD2D__ml2__p2" "petsc_options_gmres_default")
#add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/NURBS_Airfoil/TEST_Euler_NURBSAirfoil_ParametricQUAD2D__ml2__p2" "petsc_options_gmres_default")
#add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/NURBS_Airfoil/TEST_Euler_NURBSAirfoil_ParametricTRI2D__ml2__p2" "petsc_options_gmres_default")
#add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/NURBS_Airfoil/TEST_Euler_NURBSAirfoil_ParametricQUAD2D__ml2__p2" "petsc_options_gmres_default")


set (EXEC test_integration_optimization)
set (LIBS_DEPEND ${LIBS_BASE} Core Simulation Test_Integration)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
#add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/NURBS_Airfoil/TEST_Euler_NURBSAirfoil_TargetCLTestsReference_ParametricQUAD2D__ml2__p2" "petsc_options_gmres_default")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/NURBS_Airfoil/Const_DOF_NURBS/TEST_Euler_NURBSAirfoil_TargetCLTestsReference__ml1__p3" "petsc_options_gmres_default")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/NURBS_Airfoil/Const_DOF_NURBS/TEST_Euler_NURBSAirfoil_TargetCLTestsReference__ml2__p2" "petsc_options_gmres_default")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/NURBS_Airfoil/Const_DOF_NURBS/TEST_Euler_NURBSAirfoil_TargetCLTestsReference__ml3__p1" "petsc_options_gmres_default")

#- ARCHIVED OPTIMIZATION TESTS:
#add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/NURBS_Airfoil/TEST_Euler_NURBSAirfoil_ParametricQUAD2D__ml2__p1" "petsc_options_gmres_default")
#add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/NURBS_Airfoil/TEST_Euler_NURBSAirfoil_ParametricQUAD2D__ml2__p2" "petsc_options_gmres_default")
#add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/NURBS_Airfoil/TEST_Euler_NURBSAirfoil_ParametricQUAD2D__ml2__p3" "petsc_options_gmres_default")
#add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/NURBS_Airfoil/TEST_Euler_NURBSAirfoil_TargetCL_ParametricQUAD2D__ml2__p2" "petsc_options_gmres_default")
#add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/NURBS_Airfoil/TEST_Euler_NURBSAirfoil_ParametricTRI2D__ml2__p2" "petsc_options_gmres_default")
#set_tests_properties(test_integration_optimization___euler/NURBS_Airfoil/TEST_Euler_NURBSAirfoil_ParametricQUAD2D__ml2__p1 PROPERTIES TIMEOUT 7200) 
#set_tests_properties(test_integration_optimization___euler/NURBS_Airfoil/TEST_Euler_NURBSAirfoil_ParametricQUAD2D__ml2__p2 PROPERTIES TIMEOUT 13800) 
#set_tests_properties(test_integration_optimization___euler/NURBS_Airfoil/TEST_Euler_NURBSAirfoil_ParametricQUAD2D__ml2__p3 PROPERTIES TIMEOUT 13800) 

# Add tests:
# - equivalence operators.
