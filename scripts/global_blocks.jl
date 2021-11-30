using Revise

using Oiler

using CSV
using JSON
using DataFrames:eachrow
using DataFrames, DataFramesMeta
using ArchGDAL
using Setfield
const AG = ArchGDAL

using PyPlot

# options
geol_slip_rate_weight = 2.
tri_distance_weight = 20.
save_results = false


# load data
# CEA
cea_fault_file = "/home/itchy/research/gem/c_asia/c_asia_blocks/block_data/c_asia_faults.geojson"
cea_block_file = "/home/itchy/research/gem/c_asia/c_asia_blocks/block_data/c_asia_blocks.geojson"
cea_slip_rate_file = "/home/itchy/research/gem/c_asia/c_asia_blocks/block_data/c_asia_geol_slip_rates.geojson"
cea_tris_file = "/home/itchy/research/gem/c_asia/c_asia_blocks/block_data/c_asia_sub_tris.geojson"

# CHN
chn_block_file = "/home/itchy/research/gem/fault_data/china/block_data/chn_blocks.geojson"
chn_fault_file = "/home/itchy/research/gem/fault_data/china/block_data/chn_faults.geojson"
chn_slip_rate_file = "/home/itchy/research/gem/fault_data/china/block_data/geol_slip_rate_pts.geojson"

# ANA
ana_block_file = "/home/itchy/research/geodesy/global_block_comps/anatolia/block_data/anatolia_blocks.geojson"
ana_fault_file = "/home/itchy/research/geodesy/global_block_comps/anatolia/block_data/anatolia_faults.geojson"

# NEA
nea_block_file = "/home/itchy/research/geodesy/global_block_comps/ne_asia_blocks/ne_asia_blocks_.geojson"
nea_fault_file = "/home/itchy/research/geodesy/global_block_comps/ne_asia_blocks/ne_asia_faults_.geojson"
nea_slip_rate_file = "/home/itchy/research/geodesy/global_block_comps/ne_asia_blocks/ne_asia_slip_rates.geojson"

# CAS
cas_block_file = "/home/itchy/research/cascadia/cascadia_blocks/data/cascadia_blocks.geojson"
cas_fault_file = "/home/itchy/research/cascadia/cascadia_blocks/data/cascadia_block_faults.geojson"
cascadia_geol_slip_rates_file = "/home/itchy/research/cascadia/cascadia_blocks/data/cascadia_geol_slip_rate_pts.geojson"
cascadia_tris_file = "/home/itchy/research/cascadia/cascadia_blocks/data/graham_cascadia_subduction_tris.geojson"
aleut_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/alu_tris.geojson"
jdf_point_file = "/home/itchy/research/cascadia/cascadia_blocks/data/jdf_vel_pts.csv"
exp_point_file = "/home/itchy/research/cascadia/cascadia_blocks/data/explorer_vel_pts.csv"

# SUS
sus_block_file = "/home/itchy/research/us_faults/s_us_faults/s_us_blocks.geojson"
sus_fault_file = "/home/itchy/research/us_faults/s_us_faults/s_us_faults.geojson"
sus_geol_rates_file = "/home/itchy/research/us_faults/s_us_faults/new_us_faults_geol_slip_rates.geojson"

# other
glo_block_file = "/home/itchy/research/geodesy/global_block_comps/global_scale_plates/global_scale_plates.geojson"

# Geod
c_asia_gsrm_vels_file = "/home/itchy/research/gem/c_asia/c_asia_blocks/gnss_data/gsrm_c_asia_vels.geojson"
comet_gnss_vels_file = "/home/itchy/research/gem/c_asia/c_asia_blocks/gnss_data/comet_c_asia_gnss_vels.geojson"
tibet_vel_field_file = "/home/itchy/research/gem/fault_data/china/geod/tibet_vel_field.geojson"
gsrm_midas_ak_vels_file = "/home/itchy/research/cascadia/cascadia_blocks/data/vels_consolidated.geojson"



@info "joining blocks"
cea_blocks = Oiler.IO.gis_vec_file_to_df(cea_block_file)
chn_blocks = Oiler.IO.gis_vec_file_to_df(chn_block_file)
ana_blocks = Oiler.IO.gis_vec_file_to_df(ana_block_file)
nea_blocks = Oiler.IO.gis_vec_file_to_df(nea_block_file)
cas_blocks = Oiler.IO.gis_vec_file_to_df(cas_block_file)
sus_blocks = Oiler.IO.gis_vec_file_to_df(sus_block_file)
glo_blocks = Oiler.IO.gis_vec_file_to_df(glo_block_file)

block_df = vcat(cea_blocks, 
                chn_blocks,
                ana_blocks,
                nea_blocks,
                cas_blocks,
                sus_blocks,
                glo_blocks; 
                cols=:union)

println("n blocks: ", size(block_df, 1))



@info "doing faults"
asia_fault_df, asia_faults, asia_fault_vels = Oiler.IO.process_faults_from_gis_files(
                                        cea_fault_file,
                                        chn_fault_file,
                                        ana_fault_file,
                                        nea_fault_file,
                                        #cas_fault_file,
                                        #sus_fault_file;
                                        block_df=block_df)

nam_fault_df, nam_faults, nam_fault_vels = Oiler.IO.process_faults_from_gis_files(
                                        #cea_fault_file,
                                        #chn_fault_file,
                                        cas_fault_file,
                                        sus_fault_file;
                                        block_df=block_df,
                                        usd=:upper_seis_depth,
                                        lsd=:lower_seis_depth)
rename!(nam_fault_df, :upper_seis_depth => :usd)
rename!(nam_fault_df, :lower_seis_depth => :lsd)

fault_df = vcat(asia_fault_df, nam_fault_df, cols=:union)
faults = vcat(asia_faults, nam_faults)

# not sure if I need this
jdf_ridge_vels = filter( x -> x.mov == "c006", nam_fault_vels)
nam_fault_vels = filter( x -> x.mov != "c006", nam_fault_vels)

fault_vels = vcat([jdf_ridge_vels[1]], nam_fault_vels, asia_fault_vels)

println("n faults: ", length(faults))
println("n fault vels: ", length(fault_vels))


@info "doing geol slip rates"
cea_slip_rate_df = Oiler.IO.gis_vec_file_to_df(cea_slip_rate_file)
chn_slip_rate_df = Oiler.IO.gis_vec_file_to_df(chn_slip_rate_file)
nea_slip_rate_df = Oiler.IO.gis_vec_file_to_df(nea_slip_rate_file)

asia_slip_rate_df = vcat(cea_slip_rate_df, 
                         chn_slip_rate_df, 
                         nea_slip_rate_df
                         )
asia_slip_rate_df, asia_slip_rate_vels = Oiler.IO.make_geol_slip_rate_vels!(
                                            asia_slip_rate_df, fault_df)

cas_geol_slip_rate_df = Oiler.IO.gis_vec_file_to_df(cascadia_geol_slip_rates_file)
sus_geol_slip_rate_df = Oiler.IO.gis_vec_file_to_df(sus_geol_rates_file)


nam_geol_slip_rate_df = vcat(cas_geol_slip_rate_df, sus_geol_slip_rate_df)
nam_geol_slip_rate_df, nam_geol_slip_rate_vels = Oiler.IO.make_geol_slip_rate_vels!(
                                                        nam_geol_slip_rate_df, 
                                                        fault_df;
                                                        #usd="upper_seis_depth",
                                                        #lsd="lower_seis_depth",
                                                        )

geol_slip_rate_df = vcat(asia_slip_rate_df, nam_geol_slip_rate_df, cols=:union)
geol_slip_rate_vels = vcat(asia_slip_rate_vels, nam_geol_slip_rate_vels)

geol_slip_rate_df[:, :dextral_err] /= geol_slip_rate_weight
geol_slip_rate_df[:, :extension_err] /= geol_slip_rate_weight
println("n geol slip rates: ", length(geol_slip_rate_vels))




@info "doing GNSS"
cea_vel_df = Oiler.IO.gis_vec_file_to_df(c_asia_gsrm_vels_file)
com_vel_df = Oiler.IO.gis_vec_file_to_df(comet_gnss_vels_file)
tib_vel_df = Oiler.IO.gis_vec_file_to_df(tibet_vel_field_file)
cas_vel_df = Oiler.IO.gis_vec_file_to_df(gsrm_midas_ak_vels_file)

tib_vel_df[!,"station"] = string.(tib_vel_df[!,:id])


@info " doing comet gnss vels"
@time com_vels = Oiler.IO.make_vels_from_gnss_and_blocks(com_vel_df, block_df;
    ve=:v_east, vn=:v_north, ee=:sig_east, en=:sig_north, name=:name,
    fix="1111", epsg=2991,
)

@info " doing Asia gsrm vels"
@time cea_vels = Oiler.IO.make_vels_from_gnss_and_blocks(cea_vel_df, block_df;
    ve=:e_vel, vn=:n_vel, ee=:e_err, en=:n_err, name=:station,
    fix="1111", epsg=2991,
)

@info " doing comet insar vels"
tib_vels = Oiler.IO.make_vels_from_gnss_and_blocks(tib_vel_df, block_df;
    fix="1111", epsg=2991)

@info " doing NAM-rel vels"
@time cas_vels = Oiler.IO.make_vels_from_gnss_and_blocks(cas_vel_df, block_df;
                                                           fix="na",
                                                           ve=:e_vel,
                                                           vn=:n_vel,
                                                           ee=:e_err,
                                                           en=:n_err,
                                                           name=:station,
                                                           epsg=2991)

gnss_vels = vcat(com_vels, cea_vels, tib_vels, cas_vels)
println("n gnss vels: ", length(gnss_vels))


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
tris = vcat(cas_tris, mak_tris)


vels = vcat(fault_vels, 
            gnss_vels, 
            geol_slip_rate_vels, 
            jdf_vels, 
            exp_vels,
           )

vel_groups = Oiler.group_vels_by_fix_mov(vels)


@info "Solving"
@time results = Oiler.solve_block_invs_from_vel_groups(vel_groups,
            tris=tris,
            faults=faults,
            elastic_floor=1e-4,
            tri_distance_weight=tri_distance_weight,
            regularize_tris=true,
            predict_vels=true,
            pred_se=false,
            check_closures=true,
            constraint_method="kkt_sym")

Oiler.ResultsAnalysis.compare_data_results(results=results,
                                           vel_groups=vel_groups,
                                           geol_slip_rate_df=geol_slip_rate_df,
                                           geol_slip_rate_vels=geol_slip_rate_vels,
                                           fault_df=fault_df)


Oiler.WebViewer.write_web_viewer(results=results, block_df=block_df,
                                 directory="../web_viewer", ref_pole="na")


map_fig = Oiler.Plots.plot_results_map(results, vel_groups, faults, tris)
rates_fig = Oiler.Plots.plot_slip_rate_fig(geol_slip_rate_df, 
                                           geol_slip_rate_vels, 
                                           fault_df, results)

show()

println("done!")

