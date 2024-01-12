using Revise

using Oiler

using CSV
using JSON
using DataFrames:eachrow
using DataFrames, DataFramesMeta
using Setfield

using PyPlot

# options
geol_slip_rate_weight = 2.
save_results = true
pred_se = false


# load data
# CEA
cea_fault_file = "/home/itchy/research/geodesy/global_block_comps/c_asia_blocks/block_data/c_asia_faults.geojson"
cea_block_file = "/home/itchy/research/geodesy/global_block_comps/c_asia_blocks/block_data/c_asia_blocks.geojson"
cea_slip_rate_file = "/home/itchy/research/geodesy/global_block_comps/c_asia_blocks/block_data/c_asia_geol_slip_rates.geojson"
cea_tris_file = "/home/itchy/research/geodesy/global_block_comps/c_asia_blocks/block_data/c_asia_sub_tris.geojson"

# CHN
chn_block_file = "/home/itchy/research/geodesy/global_block_comps/china/block_data/chn_blocks.geojson"
chn_fault_file = "/home/itchy/research/geodesy/global_block_comps/china/block_data/chn_faults.geojson"
chn_slip_rate_file = "/home/itchy/research/geodesy/global_block_comps/china/block_data/geol_slip_rate_pts.geojson"

# ANA
ana_block_file = "/home/itchy/research/geodesy/global_block_comps/anatolia/block_data/anatolia_blocks.geojson"
ana_fault_file = "/home/itchy/research/geodesy/global_block_comps/anatolia/block_data/anatolia_faults.geojson"

# NEA
nea_block_file = "/home/itchy/research/geodesy/global_block_comps/ne_asia_blocks/ne_asia_blocks.geojson"
nea_fault_file = "/home/itchy/research/geodesy/global_block_comps/ne_asia_blocks/ne_asia_faults.geojson"
nea_slip_rate_file = "/home/itchy/research/geodesy/global_block_comps/ne_asia_blocks/ne_asia_slip_rates.geojson"
kur_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/kur_tris_slab2.geojson"

# CAS
cas_block_file = "/home/itchy/research/geodesy/global_block_comps/cascadia_blocks/data/cascadia_blocks.geojson"
cas_fault_file = "/home/itchy/research/geodesy/global_block_comps/cascadia_blocks/data/cascadia_block_faults.geojson"
cascadia_geol_slip_rates_file = "/home/itchy/research/geodesy/global_block_comps/cascadia_blocks/data/cascadia_geol_slip_rate_pts.geojson"
#cascadia_tris_file = "/home/itchy/research/cascadia/cascadia_blocks/data/graham_cascadia_subduction_tris.geojson"
cascadia_tris_file = "/home/itchy/research/geodesy/global_block_comps/cascadia_blocks/data/jdf_explorer_interface.geojson"
aleut_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/alu_tris_slab2.geojson"
#aleut_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/alu_tris.geojson"
jdf_point_file = "/home/itchy/research/geodesy/global_block_comps/cascadia_blocks/data/jdf_vel_pts.csv"
exp_point_file = "/home/itchy/research/geodesy/global_block_comps/cascadia_blocks/data/explorer_vel_pts.csv"

# SUS
sus_block_file = "/home/itchy/research/geodesy/global_block_comps/s_us_faults/s_us_blocks.geojson"
sus_fault_file = "/home/itchy/research/geodesy/global_block_comps/s_us_faults/s_us_faults.geojson"
sus_geol_rates_file = "/home/itchy/research/geodesy/global_block_comps/s_us_faults/new_us_faults_geol_slip_rates.geojson"
cali_geol_slip_rates_file = "/home/itchy/research/geodesy/global_block_comps/s_us_faults/ca_geol_slip_rates.geojson"

# SAM
sam_block_file ="/home/itchy/research/geodesy/global_block_comps/sam_blocks/block_data/sam_blocks.geojson"
sam_fault_file ="/home/itchy/research/geodesy/global_block_comps/sam_blocks/block_data/sam_faults.geojson"
sam_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/sam_slab2_60.geojson"
mora_vels_file = "/home/itchy/research/geodesy/global_block_comps/sam_blocks/block_data/mora_vels.geojson"

# CCA
cca_block_file = "/home/itchy/research/geodesy/global_block_comps/cca_blocks/block_data/cca_blocks.geojson"
cca_fault_file = "/home/itchy/research/geodesy/global_block_comps/cca_blocks/block_data/cca_faults.geojson"
cca_slip_rates_file = "/home/itchy/research/geodesy/global_block_comps/cca_blocks/block_data/cca_geol_slip_rates.geojson"
ant_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/ant_slab2.geojson"
cam_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/cam_slab2_fine.geojson"
garnier_vels_file = "/home/itchy/research/geodesy/global_block_comps/cca_blocks/geod_data/garnier_et_al_2022_vels_igs08.geojson"

# JPN
jpn_block_file = "/home/itchy/research/geodesy/global_block_comps/japan_blocks/block_data/japan_blocks.geojson"
jpn_fault_file = "/home/itchy/research/geodesy/global_block_comps/japan_blocks/block_data/japan_faults.geojson"
izu_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/izu_slab2.geojson"
ryu_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/ryu_slab2.geojson"
kjp_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/kur_jpn_slab2.geojson"

# PHL
phl_block_file = "/home/itchy/research/geodesy/global_block_comps/phil_blocks/block_data/phl_blocks.geojson"
phl_fault_file = "/home/itchy/research/geodesy/global_block_comps/phil_blocks/block_data/phl_faults.geojson"

# OCN
ocn_block_file = "/home/itchy/research/geodesy/global_block_comps/oceania_blocks/block_data/oceania_blocks.geojson"
ocn_fault_file = "/home/itchy/research/geodesy/global_block_comps/oceania_blocks/block_data/oceania_faults.geojson"
ocn_slip_rates_file = "/home/itchy/research/geodesy/global_block_comps/oceania_blocks/block_data/oceania_geol_slip_rates.geojson"
# ocn tris
png_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/png_slab2.geojson"
phi_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/phi_slab2.geojson"
man_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/man_slab2.geojson"
cot_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/cot_slab2.geojson"
sul_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/sul_slab2.geojson"
sol_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/sol_slab2.geojson"
sum_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/sum_slab2.geojson"
van_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/van_slab2.geojson"

# other
glo_block_file = "/home/itchy/research/geodesy/global_block_comps/global_scale_plates/global_scale_plates.geojson"
glo_fault_file = "/home/itchy/research/geodesy/global_block_comps/global_scale_plates/global_scale_faults.geojson"
glo_slip_rates_file = "/home/itchy/research/geodesy/global_block_comps/global_scale_plates/global_scale_slip_rates.geojson"

# Geod
c_asia_gsrm_vels_file = "/home/itchy/research/geodesy/global_block_comps/c_asia_blocks/gnss_data/gsrm_c_asia_vels.geojson"
comet_gnss_vels_file = "/home/itchy/research/geodesy/global_block_comps/c_asia_blocks/gnss_data/comet_c_asia_gnss_vels.geojson"
tibet_vel_field_file = "/home/itchy/research/geodesy/global_block_comps/china/geod/tibet_vel_field_2021_12_06.geojson"
gsrm_midas_ak_vels_file = "/home/itchy/research/geodesy/global_block_comps/cascadia_blocks/data/vels_consolidated.geojson"
weiss_vel_field_file = "/home/itchy/research/geodesy/global_block_comps/anatolia/geod_data/weiss_et_al_2020_vels_down_100.geojson"

# kur test bounds
kur_test_bounds_file = "../block_data/kur_test_bounds.geojson"
ocn_bounds_file = "../block_data/oceania_bounds.geojson"

@info "joining blocks"
cea_blocks = Oiler.IO.gis_vec_file_to_df(cea_block_file)
chn_blocks = Oiler.IO.gis_vec_file_to_df(chn_block_file)
ana_blocks = Oiler.IO.gis_vec_file_to_df(ana_block_file)
nea_blocks = Oiler.IO.gis_vec_file_to_df(nea_block_file)
cas_blocks = Oiler.IO.gis_vec_file_to_df(cas_block_file)
cca_blocks = Oiler.IO.gis_vec_file_to_df(cca_block_file)
sam_blocks = Oiler.IO.gis_vec_file_to_df(sam_block_file)
sus_blocks = Oiler.IO.gis_vec_file_to_df(sus_block_file)
jpn_blocks = Oiler.IO.gis_vec_file_to_df(jpn_block_file)
ocn_blocks = Oiler.IO.gis_vec_file_to_df(ocn_block_file)
phl_blocks = Oiler.IO.gis_vec_file_to_df(phl_block_file)
glo_blocks = Oiler.IO.gis_vec_file_to_df(glo_block_file)

block_df = vcat(cea_blocks, 
                chn_blocks,
                ana_blocks,
                nea_blocks,
                cas_blocks,
                sus_blocks,
                sam_blocks,
                cca_blocks,
                jpn_blocks,
                ocn_blocks,
                phl_blocks,
                glo_blocks; 
                cols=:union)


block_df[!, :fid] = string.(block_df[!, :fid])

println("n blocks: ", size(block_df, 1))

@info "removing antarctica for now"
ant_df = filter(row -> (row.fid == "ant"), block_df)
block_df = filter(row -> !(row.fid == "ant"), block_df)

# turn this on to precompile
#bound_df = Oiler.IO.gis_vec_file_to_df(kur_test_bounds_file)
#block_df = Oiler.IO.get_blocks_in_bounds!(block_df, bound_df; epsg=102016)
#println("n blocks: ", size(block_df, 1))

@info "doing GNSS"
cea_vel_df = Oiler.IO.gis_vec_file_to_df(c_asia_gsrm_vels_file)
com_vel_df = Oiler.IO.gis_vec_file_to_df(comet_gnss_vels_file)
tib_vel_df = Oiler.IO.gis_vec_file_to_df(tibet_vel_field_file)
cas_vel_df = Oiler.IO.gis_vec_file_to_df(gsrm_midas_ak_vels_file)
gar_vel_df = Oiler.IO.gis_vec_file_to_df(garnier_vels_file)
mora_vel_df = Oiler.IO.gis_vec_file_to_df(mora_vels_file)
weiss_vel_field_df = Oiler.IO.gis_vec_file_to_df(weiss_vel_field_file)

tib_vel_df[!,"station"] = string.(tib_vel_df[!,:fid])
weiss_vel_field_df[!,"station"] = map(x->join(["weiss_", x]), 
                                      string.(weiss_vel_field_df[!,:fid]))


@info " doing comet gnss vels"
@time com_vels = Oiler.IO.make_vels_from_gnss_and_blocks(com_vel_df, block_df;
    ve=:v_east, vn=:v_north, ee=:sig_east, en=:sig_north, name=:name,
    fix="1111", epsg=102016,
)

@info " doing Asia gsrm vels"
@time cea_vels = Oiler.IO.make_vels_from_gnss_and_blocks(cea_vel_df, block_df;
    ve=:e_vel, vn=:n_vel, ee=:e_err, en=:n_err, name=:station,
    fix="1111", epsg=102016,
)

@info " doing comet insar vels"
@time tib_vels = Oiler.IO.make_vels_from_gnss_and_blocks(tib_vel_df, block_df;
    fix="1111", epsg=102016)

@info " doing weiss insar vels"
@time wss_vels = Oiler.IO.make_vels_from_gnss_and_blocks(
    weiss_vel_field_df, block_df;
    ve=:e_vel, vn=:n_vel, ee=:e_err, en=:n_err, name=:station,
    fix="1111"
)

@info " doing Garnier vels"
@time gar_vels = Oiler.IO.make_vels_from_gnss_and_blocks(gar_vel_df, block_df;
    fix="igs08", epsg=102016,
    ve=:e_vel, vn=:n_vel, ee=:e_rr, en=:n_err, name=:site
)

@time mora_vels = Oiler.IO.make_vels_from_gnss_and_blocks(mora_vel_df, block_df;
    fix="itrf14", epsg=102016,
    ve=:e_vel, vn=:n_vel, ee=:e_err, en=:n_err, name=:station
)
@info " doing NAM-rel vels"
@time cas_vels = Oiler.IO.make_vels_from_gnss_and_blocks(cas_vel_df, block_df;
                                                           fix="na",
                                                           ve=:e_vel,
                                                           vn=:n_vel,
                                                           ee=:e_err,
                                                           en=:n_err,
                                                           name=:station,
                                                           epsg=102016)


@info " doing NAM-rel vels (Antarctica)"
@time ant_vels = Oiler.IO.make_vels_from_gnss_and_blocks(cas_vel_df, ant_df;
                                                           fix="1111",
                                                           ve=:e_vel,
                                                           vn=:n_vel,
                                                           ee=:e_err,
                                                           en=:n_err,
                                                           name=:station,
                                                           epsg=102019)

@info "re-adding Antarctica to blocks"
block_df = vcat(block_df, ant_df)

gnss_vels = vcat(com_vels, 
                 cea_vels, 
                 tib_vels,
                 gar_vels,
                 cas_vels,
                 ant_vels,
                 wss_vels,
                 )

println("n gnss vels: ", length(gnss_vels))
@info "doing faults"
asia_fault_df, asia_faults, asia_fault_vels = Oiler.IO.process_faults_from_gis_files(
                                        cea_fault_file,
                                        chn_fault_file,
                                        ana_fault_file,
                                        nea_fault_file,
                                        sam_fault_file,
                                        cca_fault_file,
                                        glo_fault_file,
                                        ocn_fault_file,
                                        phl_fault_file,
                                        block_df=block_df,
                                        subset_in_bounds=true,
                                        check_blocks=false)


# need to have large uncertaintes for Japan
jpn_fault_df, jpn_faults, jpn_fault_vels = Oiler.IO.process_faults_from_gis_files(
                                        jpn_fault_file,
                                        block_df=block_df,
                                        usd_default=1.,
                                        lsd_default=10.,
                                        e_default=1.,
                                        check_blocks=false,
                                        )


nam_fault_df, nam_faults, nam_fault_vels = Oiler.IO.process_faults_from_gis_files(
                                        cas_fault_file,
                                        sus_fault_file;
                                        block_df=block_df,
                                        usd=:upper_seis_depth,
                                        lsd=:lower_seis_depth,
                                        subset_in_bounds=true,
                                        check_blocks=false,
                                        )
rename!(nam_fault_df, :upper_seis_depth => :usd)
rename!(nam_fault_df, :lower_seis_depth => :lsd)

fault_df = vcat(asia_fault_df,
                nam_fault_df, 
                jpn_fault_df, 
                cols=:union)
faults = vcat(asia_faults, 
              nam_faults, 
              jpn_faults,
              )

# not sure if I need this
#jdf_ridge_vels = filter( x -> x.mov == "c006", nam_fault_vels)
#nam_fault_vels = filter( x -> x.mov != "c006", nam_fault_vels)

fault_vels = vcat(#[jdf_ridge_vels[1]], 
                  nam_fault_vels, 
                  asia_fault_vels,
                  jpn_fault_vels,
                  )

println("n faults: ", length(faults))
println("n fault vels: ", length(fault_vels))


@info "checking fault hw and fw"
ffs = Oiler.Utils.check_hw_fw_all(faults, block_df; verbose=true)

println("n faults left: ", length(ffs))


@info "doing non-fault block boundaries"
@time non_fault_bounds = Oiler.IO.get_non_fault_block_bounds(block_df, faults, verbose=true)
bound_vels = vcat(map(x->Oiler.Boundaries.boundary_to_vels(x, ee=3.0, en=3.0), non_fault_bounds)...)
println("n non-fault-bound vels: ", length(bound_vels))



@info "doing geol slip rates"
cea_slip_rate_df = Oiler.IO.gis_vec_file_to_df(cea_slip_rate_file)
chn_slip_rate_df = Oiler.IO.gis_vec_file_to_df(chn_slip_rate_file)
nea_slip_rate_df = Oiler.IO.gis_vec_file_to_df(nea_slip_rate_file)
cca_slip_rate_df = Oiler.IO.gis_vec_file_to_df(cca_slip_rates_file)
ocn_slip_rate_df = Oiler.IO.gis_vec_file_to_df(ocn_slip_rates_file)
glo_slip_rate_df = Oiler.IO.gis_vec_file_to_df(glo_slip_rates_file)

asia_slip_rate_df = vcat(cea_slip_rate_df, 
                         chn_slip_rate_df, 
                         nea_slip_rate_df,
                         cca_slip_rate_df,
                         ocn_slip_rate_df,
                         glo_slip_rate_df,
                         )
asia_slip_rate_df, asia_slip_rate_vels = Oiler.IO.make_geol_slip_rate_vels!(
                                            asia_slip_rate_df, fault_df)

cas_geol_slip_rate_df = Oiler.IO.gis_vec_file_to_df(cascadia_geol_slip_rates_file)
sus_geol_slip_rate_df = Oiler.IO.gis_vec_file_to_df(sus_geol_rates_file)
cal_geol_slip_rate_df = Oiler.IO.gis_vec_file_to_df(cali_geol_slip_rates_file)


nam_geol_slip_rate_df = vcat(cas_geol_slip_rate_df, sus_geol_slip_rate_df, cal_geol_slip_rate_df)
nam_geol_slip_rate_df, nam_geol_slip_rate_vels = Oiler.IO.make_geol_slip_rate_vels!(
                                                        nam_geol_slip_rate_df, 
                                                        fault_df;
                                                        weight=geol_slip_rate_weight,
                                                        #usd="upper_seis_depth",
                                                        #lsd="lower_seis_depth",
                                                        )

geol_slip_rate_df = vcat(asia_slip_rate_df, 
                         nam_geol_slip_rate_df, 
                         cols=:union)
geol_slip_rate_vels = vcat(asia_slip_rate_vels, 
                           nam_geol_slip_rate_vels)


println("n geol slip rates: ", length(geol_slip_rate_vels))






@info "doing Explorer and JdF plates"

# Fake JDF vel points
jdf_pts = CSV.read(jdf_point_file, DataFrame)

jdf_na_pole = Oiler.PoleSphere(
    lon=68.3,
    lat=-32.0,
    rotrate=1.081,
    elat=0.5,
    elon=0.5,
    erotrate=0.15,
    fix="na",
    mov="c006")


jdf_vels = Oiler.BlockRotations.predict_block_vels(jdf_pts[:,:lon], 
    jdf_pts[:,:lat], jdf_na_pole)


jdf_vels = [Oiler.VelocityVectorSphere(vel; vel_type="fault") for vel in jdf_vels]
#jdf_vels = jdf_vels


exp_pac_pole = Oiler.PoleSphere(
    lat=54.80, lon=-116.62, rotrate=3.1, 
    elat=0.25, elon=0.25, erotrate=0.31,
    fix="c024", mov="c112")

exp_pts = CSV.read(exp_point_file, DataFrame)

exp_vels = Oiler.BlockRotations.predict_block_vels(exp_pts[:,:lon],
                                                   exp_pts[:,:lat],
                                                   exp_pac_pole)

exp_vels = [Oiler.VelocityVectorSphere(vel; vel_type="fault") for vel in exp_vels]

@info "doing tris"
cas_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(cascadia_tris_file))
mak_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(cea_tris_file))
alu_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(aleut_tris_file))
kur_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(kur_tris_file))
kjp_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(kjp_tris_file))
izu_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(izu_tris_file))
ryu_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(ryu_tris_file))


pac_eur_pole = Oiler.PoleCart(
  x = -1.1959994895724632e-9,
  y = 7.658839979795385e-9,
  z = -1.4175151626090377e-8,
  ex = 2.9245747805725353e-11,
  ey = 1.8371290128703008e-11,
  ez = 1.8786487072346268e-11,
  cxy = 2.9286193076352114e-22,
  cxz = -1.6021810026353652e-22,
  cyz = -9.176905987901183e-23,
  fix = "1111",
  mov = "c024"
)

@info "kur tris"
kur_tris = Oiler.Utils.tri_priors_from_pole(kur_tris, pac_eur_pole,
                                            locking_fraction=0.8,
                                            depth_adjust=true,
                                            err_coeff=10000.)

@info "kjp tris"
kjp_tris = Oiler.Utils.tri_priors_from_pole(kjp_tris, pac_eur_pole,
                                            locking_fraction=1.,
                                            depth_adjust=true,
                                            err_coeff=10000.)



pac_phi_pole = Oiler.PoleCart(
  x = -1.3310752536689726e-8,
  y = 1.310215442152588e-8,
  z = 1.9121268659160096e-9,
  ex = 1.2623134821110894e-10,
  ey = 1.1309767520167837e-10,
  ez = 9.90699121703064e-11,
  cxy = -1.2383643596201617e-20,
  cxz = -1.1007058896261733e-20,
  cyz = 9.643144806376538e-21,
  fix = "phi",
  mov = "pac"
)

@info "izu tris"
izu_tris = Oiler.Utils.tri_priors_from_pole(izu_tris, pac_phi_pole,
                                            locking_fraction=0.8,
                                            depth_adjust=true,
                                            err_coeff=10000.)

phi_eur_pole = Oiler.Utils.PoleCart(
  x = 1.2114753047117263e-8,
  y = -5.4433144417304956e-9,
  z = -1.6087278492006386e-8,
  ex = 1.227967406183004e-10,
  ey = 1.1159560849348577e-10,
  ez = 9.727237737874402e-11,
  cxy = -1.2676505526965139e-20,
  cxz = -1.0846840795998196e-20,
  cyz = 9.734913866255549e-21,
  fix = "1111",
  mov = "phi",
)

@info "ryu tris"
ryu_tris = Oiler.Utils.tri_priors_from_pole(ryu_tris, phi_eur_pole,
                                            locking_fraction=0.8,
                                            depth_adjust=true,
                                            err_coeff=10000.)


pac_na_pole = Oiler.PoleCart(
  x = -2.4885570284209037e-9,
  y = 8.602351086527437e-9,
  z = -1.0424548581912783e-8,
  ex = 5.4017557235402144e-11,
  ey = 3.829966190320154e-11,
  ez = 4.100348456652813e-11,
  cxy = 1.0584250692981392e-21,
  cxz = -9.478742893924508e-22,
  cyz = -7.463051012669244e-22,
  fix = "na",
  mov = "c024",
)

@info "alu tris"
alu_tris = Oiler.Utils.tri_priors_from_pole(alu_tris, pac_na_pole,
                                              locking_fraction=0.5,
                                              depth_adjust=true,
                                              err_coeff=10000.)


function set_tri_rates(tri; ds=20., de=10000., ss=0., se=10000., frac=0.5)
    tri = @set tri.dip_slip_rate = ds * frac
    tri = @set tri.dip_slip_err = de * frac
    tri = @set tri.strike_slip_rate = ss
    tri = @set tri.strike_slip_err = se
    tri
end

@info "mak tris"
mak_tris = map(set_tri_rates, mak_tris)

@info "ant tris"
ant_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(ant_tris_file))
ant_tris = map(x->set_tri_rates(x; ds=5.), ant_tris)

@info "cam tris"
cam_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(cam_tris_file))
cam_tris = map(set_tri_rates, cam_tris)

@info "sam tris"
sam_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(sam_tris_file))
sam_tris = map(set_tri_rates, sam_tris)


# ocn tris, priors from Bird 2002 (roughly), 30% coupling
sol_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(sol_tris_file))
sul_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(sul_tris_file))
sum_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(sum_tris_file))
png_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(png_tris_file))
van_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(van_tris_file))

@info "sol tris"
sol_tris = map(x->set_tri_rates(x; ds=85., ss=-40., frac=0.3), sol_tris)
sul_tris = map(x->set_tri_rates(x; ds=38., ss=-14., frac=0.3), sul_tris)
sum_tris = map(x->set_tri_rates(x; ds=35., ss=36.,  frac=0.3), sum_tris)
png_tris = map(x->set_tri_rates(x; ds=27., ss=2., frac=0.3), png_tris)
van_tris = map(x->set_tri_rates(x; ds=100., ss=-13., frac=0.3), van_tris)


tris = vcat(cas_tris, 
            mak_tris, 
            alu_tris, 
            kur_tris,
            kjp_tris,
            izu_tris,
            ryu_tris,
            ant_tris,
            cam_tris,
            sam_tris,
            sol_tris,
            sul_tris,
            sum_tris,
            png_tris,
            van_tris
           )

println("n tris: ", length(tris) )

vels = vcat(fault_vels, 
            gnss_vels, 
            geol_slip_rate_vels, 
            jdf_vels, 
            exp_vels,
            bound_vels,
           )

vel_groups = Oiler.group_vels_by_fix_mov(vels)

tri_distance_weight = 5.

@info "Solving"
@time results = Oiler.solve_block_invs_from_vel_groups(vel_groups,
            tris=tris,
            faults=faults,
            #elastic_floor=0.,
            elastic_floor=1e-4,
            tri_distance_weight=tri_distance_weight,
            regularize_tris=true,
            tri_priors=true,
            predict_vels=true,
            pred_se=pred_se,
            se_iters=200,
            check_closures=false,
            constraint_method="kkt_sym",
            check_nans=true,
            sparse_lhs=true,
            factorization="lu")

Oiler.ResultsAnalysis.compare_data_results(results=results,
                                           vel_groups=vel_groups,
                                           geol_slip_rate_df=geol_slip_rate_df,
                                           geol_slip_rate_vels=geol_slip_rate_vels,
                                           fault_df=fault_df)


if save_results == true
    Oiler.IO.write_tri_results_to_gj(tris, results,
                                     "../results/all_tris_no_errs.geojson",
                                     name="global tri results")
    
    if pred_se
        Oiler.IO.write_fault_results_to_gj(results,
                                       "../results/global_faults.geojson",
                                       name="global fault results",
                                       calc_rake=true, calc_slip_rate=true)
    else
        Oiler.IO.write_fault_results_to_gj(results,
                                       "../results/global_faults_no_errs.geojson",
                                       name="global fault results",
                                       calc_rake=true, calc_slip_rate=true)
    end
    Oiler.IO.write_gnss_vel_results_to_csv(results, vel_groups;
                                       name="../results/global_gnss_results.csv")
end

Oiler.WebViewer.write_web_viewer(results=results, block_df=block_df,
                                 directory="../web_viewer", ref_pole="na")


map_fig = Oiler.Plots.plot_results_map(results, vel_groups, faults, tris)
rates_fig = Oiler.Plots.plot_slip_rate_fig(geol_slip_rate_df, 
                                           geol_slip_rate_vels, 
                                           fault_df, results)

show()

println("done!")

