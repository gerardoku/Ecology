# Script para simulación de dinámicas de ocupación en una metapoblación con paisajes simulados usando el paquete MetaLandSim
# from Mestre et al. 2018

# La idea sería usar este paquete para escribir un paper de las dinámicas de ocupación y recolonización después de la deforestación
# También se podría hacer uso de este paquete para simular lo que pasó antes de la actualidad con la deforestación histórica ****

# install.packages("MetaLandSim")
library(MetaLandSim)

# Funciones básicas ----

par(mfrow = c(1,2))
 rl <- rland.graph(mapsize = 1000, 
                   dist_m = 60,
                   areaM = 0.5, 
                   areaSD = 0.2, 
                   Npatch =70, 
                   disp = 100, 
                   plotG = TRUE)

 sp_t0 <- species.graph(rl=rl, 
                        method="percentage", 
                        parm = 50, 
                        nsew="none", 
                        plotG=TRUE)


# Simulation with graphs (metapops) ----

# Although this procedure can be carried out using the functions mentioned above,
# it is easier to complete the full simulation using only one function that runs all
# the others internally, while allowing for a repetition of the process, the function
# is ’iterate.graph’.
# Here the simulation process will run only with 2 iterations, for demonstration:

# load species parameters 
 data(param1)

# run simulations
 it1 <- iterate.graph(iter = 2,
                     mapsize = 1000,
                     dist_m = 30,
                     areaM = 0.5,
                     areaSD= 0.1,
                     Npatch = 200,
                     disp = 800,
                     span = 100,
                     par1 = "stoc",
                     par2 = 2,
                     par3 = 2,
                     method = "percentage",
                     parm = 50,
                     nsew = "none",
                     succ = "none",
                     param_df = param1,
                     kern = "op1",
                     conn = "op1",
                     colnz = "op1",
                     ext = "op1",
                     b = 1,
                     graph = TRUE) # Si pongo TRUE se genera un archivo HTML con los gráficos de las distribuciones

# As a result the user will have a large number of simulations which represent
# the occupation of a species with a given set of characteristics (as defined by
# the parameters) in a dynamic landscape. The advantage of this approach is
# that it requires less parameters which can be estimated from real occupancy or
# turnover data. It does not require demographic data, the parameters can be
# derived using only patch occupancy data of one snapshot or sampling session
# (ideally more).
# Here the simulation procedure is repeated only twice (parameter ’iter’),
# although more simulations have to be run in order to obtain robust results.
# However, depending on computing power, this simulation can be highly timeconsuming
# (from hours to several days).
# After version 0.5 of MetaLandSim an aditional option was made available to
# the users: the argument ’succ’. This allows to chose different options regarding
# the species preference relating the successional stage of habitat patches, with the
# following options: ’none’ - No discrimination regarding patch successional stage;
# ’early’ - The species prefers patches in an earlier successional stage; ’mid’ - The
# species prefers patches in the mid of the succession; ’late’ - The species prefers
# patches in an later successional stage. This new option includes an additional
# factor to the extinction probability, changing it with patch age, as follows:
 
# figure missing

# In the above figure it is visible that, when succ=’early’, the extinction concerning
# successional stage is lower for younger patches (left); when succ=’mid’,
# the extinction concerning successional stage is lower for habitat patches with intermediate
# ages; when succ=’late’, the extinction concerning successional stage
# is lower for older patches. This extinction factor is combined with the one derived
# from the application of the SPOM model considered (generally depending
#                                                    on patch area).
