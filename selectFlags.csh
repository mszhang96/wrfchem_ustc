#!/bin/csh
 

set inFile = $1
#physics
set filesWeCare_physics = ( \
module_ra_rrtmg_lw.f90 \
module_ra_rrtmg_sw.f90 \
module_radiation_driver.f90 \
module_microphysics_driver.f90 \
module_mp_morr_two_moment.f90 \
module_physics_addtendc.f90 \
module_physics_init.f90 \
)
#chemistry
set filesWeCare_chem = (  \
aerosol_driver.f90  \
chem_driver.f90     \
chemics_init.f90    \
cloudchem_driver.f90  \
convert_emiss.f90  \
dry_dep_driver.f90  \
emissions_driver.f90 \
mechanism_driver.f90  \
module_add_emis_cptec.f90 \
module_add_emiss_burn.f90 \
module_aer_drydep.f90  \
module_aer_opt_out.f90  \
module_aerosols_sorgam.f90  \
module_bioemi_beis313.f90  \
module_bioemi_megan2.f90  \
module_bioemi_simple.f90  \
module_cbm4_addemiss.f90  \
module_cbm4_initmixrats.f90  \
module_cbmz_addemiss.f90   \
module_cbmz.f90   \
module_cbmz_initmixrats.f90   \
module_cbmz_lsodes_solver.f90   \
module_cbmz_rodas3_solver.f90   \
module_cbmz_rodas_prep.f90   \
module_chem_plumerise_scalar.f90   \
module_chem_utilities.f90   \
module_cmu_bulkaqchem.f90   \
module_cmu_dvode_solver.f90   \
module_convtrans_prep.f90   \
module_ctrans_aqchem.f90  \
module_ctrans_grell.f90  \
module_ctrans_kfeta.f90  \
module_data_cbmz.f90  \
module_data_cmu_bulkaqchem.f90    \
module_data_gocartchem.f90    \
module_data_gocart_dust.f90   \
module_data_gocart_seas.f90   \
module_data_isrpia.f90    \
module_data_megan2.f90    \
module_data_mgn2mech.f90   \
module_data_mosaic_asect.f90   \
module_data_mosaic_other.f90   \
module_data_mosaic_therm.f90   \
module_data_radm2.f90   \
module_data_rrtmgaeropt.f90   \
module_data_sorgam.f90   \
module_dep_simple.f90   \
module_emissions_anthropogenics.f90   \
module_fastj_data.f90    \
module_fastj_mie.f90    \
module_ftuv_driver.f90   \
module_ftuv_subs.f90   \
module_gocart_aerosols.f90   \
module_gocart_chem.f90   \
module_gocart_dmsemis.f90  \
module_gocart_drydep.f90  \
module_gocart_dust.f90  \
module_gocart_seasalt.f90   \
module_gocart_settling.f90  \
module_gocart_so2so4.f90   \
module_hetn2o5.f90   \
module_input_chem_bioemiss.f90    \
module_input_chem_data.f90    \
module_input_dust_errosion.f90   \
module_input_gocart_dms.f90   \
module_input_tracer_data.f90   \
module_input_tracer.f90   \
module_isrpia.f90  \
module_lightning_driver.f90   \
module_ltng_crm.f90   \
module_mixactivate_wrappers.f90  \
module_mosaic_addemiss.f90   \
module_mosaic_cloudchem.f90   \
module_mosaic_coag.f90   \
module_mosaic_csuesat.f90   \
module_mosaic_driver.f90   \
module_mosaic_drydep.f90   \
module_mosaic_initmixrats.f90   \
module_mosaic_movesect.f90  \
module_mosaic_newnuc.f90    \
module_mosaic_therm.f90    \
module_mosaic_wetscav.f90   \
module_optical_averaging.f90   \
module_peg_util.f90    \
module_phot_fastj.f90    \
module_phot_mad.f90    \
module_plumerise1.f90    \
module_radm.f90   \
module_sea_salt_emis.f90   \
module_sorgam_cloudchem.f90  \
module_vash_settling.f90   \
module_vertmx_wrf.f90   \
module_wave_data.f90   \
module_wetdep_ls.f90   \
module_wetscav_driver.f90   \
module_zero_plumegen_coms.f90   \
optical_driver.f90   \
photolysis_driver.f90            )

#Dynamics
set filesWeCare_dyn = (module_first_rk_step_part1.f90    \
        module_first_rk_step_part2.f90    \
    	solve_em.f90  )

#Gather all files
set filesWeCare = ( $filesWeCare_physics $filesWeCare_chem $filesWeCare_dyn )


set flags = '-g' #For files we dont care about
foreach file ( $filesWeCare )

    if ( $inFile == $file ) then
	 set flags = '--O0 --nfix -g  --trace --trap --chk a,s,u --varheap --pca' #For files we care about
	 break
    endif
end
echo $flags    
