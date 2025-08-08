#BASELINE SCENARIO
#Reproduce a situation with high solar irradiance and low residential demand

#Strategy: combined active and reactive power control
#An 18-node test network is used

######### 1.SYSTEM DATA ##############

nN=18 #Number of nodes

#Load data (in W and var)
P_Demand_a = [0, 1500, 0, 0, 0, 900, 1100, 0, 0, 0, 1100, 600, 0, 0, 900, 0, 0, 500]  # Active power demanded on phase a of each node
P_Demand_b = [0, 0, 1100, 0, 0, 0, 0, 600, 0, 0, 800, 0, 900, 0, 600, 0, 1000, 700]   # Active power demanded on phase b of each node
P_Demand_c = [0, 0, 0, 1000, 0, 0, 0, 0, 900, 0, 1400, 700, 0, 1100, 0, 0, 1000, 0]   # Active power demanded on phase c of each node
Q_Demand_a = [0, 130, 0, 0, 0, 50, 40, 0, 0, 0, 100, 90, 0, 0, 50, 0, 0, 90]          # Reactive power demanded on phase a of each node
Q_Demand_b = [0, 0, 50, 0, 0, 0, 0, 60, 0, 0, 40, 0, 100, 0, 80, 0, 100, 90]          # Reactive power demanded on phase b of each node
Q_Demand_c = [0, 0, 0, 20, 0, 0, 0, 0, 60, 0, 100, 150, 0, 60, 0, 0, 30, 0]           # Reactive power demanded on phase c of each node

#Conventional Generation Data (in W and var)
P_Cost_G = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Generation cost in each node
P_Gen_lb = [-300000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Lower limit of active generation of each node (for all 3 phases)
P_Gen_ub = [300000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]    # Upper limit of active generation of each node (for all 3 phases)
Q_Gen_lb = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]    # Lower limit of reactive generation of each node (for all 3 phases)
Q_Gen_ub = [50000,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Upper limit of reactive generation of each node (for all 3 phases)
Gen_Status = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]   # Generation status (1: conventional generator present at the node; 0: not present)


#Maximum available apparent power of the photovoltaic generators (in VA)
S_Sol_Gen_ub_a = [0, 2500, 0, 0, 0, 2500, 0, 0, 0, 0, 2500, 2500, 0, 0, 0, 0, 0, 2500]        # In phase a
S_Sol_Gen_ub_b = [0, 0, 2500, 0, 0, 0, 0, 2500, 0, 0, 2500, 0, 0, 0, 2500, 0, 2500, 2500]     # In phase b
S_Sol_Gen_ub_c = [0, 0, 0, 0, 0, 0, 0, 0, 2500, 0, 2500, 2500, 0, 2500, 0, 0, 2500, 0]        # In phase c

#Node voltage limits
#Phase-to-neutral voltages (in V)
# ±3% limits are applied
Vmin=0.97*400/sqrt(3)
Vmax=1.03*400/sqrt(3)
V_Nodo_lb=[Vmin, Vmin, Vmin, Vmin, Vmin, Vmin, Vmin, Vmin, Vmin, Vmin, Vmin, Vmin, Vmin, Vmin, Vmin, Vmin, Vmin, Vmin]  # Lower voltage limit (per node, common to all three phases)
V_Nodo_ub=[Vmax, Vmax, Vmax, Vmax, Vmax, Vmax, Vmax, Vmax, Vmax, Vmax, Vmax, Vmax, Vmax, Vmax, Vmax, Vmax, Vmax, Vmax]  # Upper voltage limit (per node, common to all three phases)


#Line data
nL = 17  # Number of lines
dLinea_fbus = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 5, 13, 14, 15, 16, 16]     # Starting node of the line
dLinea_tbus = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]   # End node of the line

#Hypothesis: all lines have the same parameters (R1=0.446 X1=0.071 R0=1.505 X0=0.083 C1=0 C0=0 in ohm/km) and same lenght (0.1 km)

# Series impedance matrix of a line (in ohm)
Zij=0.1 .* [0.799+im*0.075 0.353+im*0.004 0.353+im*0.004; 0.353+im*0.004 0.799+im*0.075 0.353+im*0.004; 0.353+im*0.004 0.353+im*0.004 0.799+im*0.075]   #(en ohmios)

#Shunt admittance matrix of a line (in ohm)
Ysh=[0 0 0; 0 0 0; 0 0 0] #(en ohmios)

# Admittance matrix of the entire system
# It is built by assembling the contributions of each line
N=[0 0 0; 0 0 0; 0 0 0] # 3x3 null matrix
a=inv(Zij)
Y=[ 
a -a N N N N N N N N N N N N N N N N;
-a 2*a -a N N N N N N N N N N N N N N N;
N -a 2*a -a N N N N N N N N N N N N N N;
N N -a 2*a -a N N N N N N N N N N N N N;
N N N -a 3*a -a N N N N N N -a N N N N N;
N N N N -a 2*a -a N N N N N N N N N N N;
N N N N N -a 2*a -a N N N N N N N N N N;
N N N N N N -a 2*a -a N N N N N N N N N;
N N N N N N N -a 2*a -a N N N N N N N N;
N N N N N N N N -a 3*a -a -a N N N N N N;
N N N N N N N N N -a a N N N N N N N;
N N N N N N N N N -a N a N N N N N N;
N N N N -a N N N N N N N 2*a -a N N N N;
N N N N N N N N N N N N -a 2*a -a N N N;
N N N N N N N N N N N N N -a 2*a -a N N;
N N N N N N N N N N N N N N -a 3*a -a -a;
N N N N N N N N N N N N N N N -a a N;
N N N N N N N N N N N N N N N -a N a;
]


######### 2.UBOPF MODEL AND SOLUTION ############

using JuMP, Ipopt # Package loading

# Solver configuration
model = Model(Ipopt.Optimizer)  # The problem is solved using the Ipopt solver 
set_optimizer_attribute(model, "max_iter", 5000000) # Maximum number of iterations
set_optimizer_attribute(model, "tol", 1e-7)  # Tolerance
set_silent(model)   # Disable default solver output

# VARIABLE DEFINITION

#Active and reactive power of conventional generators (per phase)
@variable(model, P_G_a[i in 1:nN])
@variable(model, Q_G_a[i in 1:nN])
@variable(model, P_G_b[i in 1:nN])
@variable(model, Q_G_b[i in 1:nN])
@variable(model, P_G_c[i in 1:nN])
@variable(model, Q_G_c[i in 1:nN])

#Active and reactive power of photovoltaic generators (per phase)
@variable(model, P_G_solar_a[i in 1:nN])
@variable(model, P_G_solar_b[i in 1:nN])
@variable(model, P_G_solar_c[i in 1:nN])
@variable(model, Q_G_solar_a[i in 1:nN])
@variable(model, Q_G_solar_b[i in 1:nN])
@variable(model, Q_G_solar_c[i in 1:nN])


#Nodal voltages in rectangular form (real and imaginary parts)
@variable(model, V_a_r[1:nN], start=(400/sqrt(3))*cos(0))
@variable(model, V_a_i[1:nN], start=(400/sqrt(3))*sin(0))
@variable(model, V_b_r[1:nN], start=(400/sqrt(3))*cos(-2pi/3))
@variable(model, V_b_i[1:nN], start=(400/sqrt(3))*sin(-2pi/3))
@variable(model, V_c_r[1:nN], start=(400/sqrt(3))*cos(2pi/3))
@variable(model, V_c_i[1:nN], start=(400/sqrt(3))*sin(2pi/3))


#Construction of the complex voltage vector per phase
v_fase = []
for i in 1:nN
    push!(v_fase, V_a_r[i] + im*V_a_i[i])
    push!(v_fase, V_b_r[i] + im*V_b_i[i])
    push!(v_fase, V_c_r[i] + im*V_c_i[i])
end


# OBJECTIVE FUNCTION
#Minimization of the total active power generation cost of conventional generators
@objective(model, Min, sum(P_Cost_G[i]*(P_G_a[i] + P_G_b[i] + P_G_c[i]) for i in 1:nN))

# CONSTRAINTS

# Voltage reference at node 1 
@constraint(model, V_a_r[1]==(400/sqrt(3))*cos(0))
@constraint(model, V_a_i[1]==(400/sqrt(3))*sin(0))
@constraint(model, V_b_r[1]==(400/sqrt(3))*cos(-2pi/3))
@constraint(model, V_b_i[1]==(400/sqrt(3))*sin(-2pi/3))
@constraint(model, V_c_r[1]==(400/sqrt(3))*cos(2pi/3))
@constraint(model, V_c_i[1]==(400/sqrt(3))*sin(2pi/3))

# Generation limits of the conventional generators
@constraint(model, [i in 1:nN], P_Gen_lb[i] * Gen_Status[i] <= P_G_a[i] <= P_Gen_ub[i] * Gen_Status[i])
@constraint(model, [i in 1:nN], Q_Gen_lb[i] * Gen_Status[i] <= Q_G_a[i] <= Q_Gen_ub[i] * Gen_Status[i])
@constraint(model, [i in 1:nN], P_Gen_lb[i] * Gen_Status[i] <= P_G_b[i] <= P_Gen_ub[i] * Gen_Status[i])
@constraint(model, [i in 1:nN], Q_Gen_lb[i] * Gen_Status[i] <= Q_G_b[i] <= Q_Gen_ub[i] * Gen_Status[i])
@constraint(model, [i in 1:nN], P_Gen_lb[i] * Gen_Status[i] <= P_G_c[i] <= P_Gen_ub[i] * Gen_Status[i])
@constraint(model, [i in 1:nN], Q_Gen_lb[i] * Gen_Status[i] <= Q_G_c[i] <= Q_Gen_ub[i] * Gen_Status[i])

# Generation limits of the photovoltaic generators
@constraint(model, [i in 1:nN], (P_G_solar_a[i])^2 + (Q_G_solar_a[i])^2 <= (S_Sol_Gen_ub_a[i])^2)
@constraint(model, [i in 1:nN], (P_G_solar_b[i])^2 + (Q_G_solar_b[i])^2 <= (S_Sol_Gen_ub_b[i])^2)
@constraint(model, [i in 1:nN], (P_G_solar_c[i])^2 + (Q_G_solar_c[i])^2 <= (S_Sol_Gen_ub_c[i])^2)

# Nodal voltage limits
@constraint(model, [i in 1:nN], (V_Nodo_lb[i])^2 <= (V_a_r[i])^2 +(V_a_i[i])^2 <= (V_Nodo_ub[i])^2)
@constraint(model, [i in 1:nN], (V_Nodo_lb[i])^2 <= (V_b_r[i])^2 +(V_b_i[i])^2 <= (V_Nodo_ub[i])^2)
@constraint(model, [i in 1:nN], (V_Nodo_lb[i])^2 <= (V_c_r[i])^2 +(V_c_i[i])^2 <= (V_Nodo_ub[i])^2)

# Active and reactive power nodal balance
@constraint(model, [i in 1:nN], [P_G_a[i]; P_G_b[i]; P_G_c[i]] + [P_G_solar_a[i]; P_G_solar_b[i]; P_G_solar_c[i]] - [P_Demand_a[i]; P_Demand_b[i]; P_Demand_c[i]] == real([V_a_r[i] + im*(V_a_i[i]); V_b_r[i] + im*(V_b_i[i]); V_c_r[i] + im*(V_c_i[i])] .* conj([conj(Y[3*i-2,:]')*v_fase; conj(Y[3*i-1,:]')*v_fase; conj(Y[3*i,:]')*v_fase]) ))
@constraint(model, [i in 1:nN], [Q_G_a[i]; Q_G_b[i]; Q_G_c[i]] + [Q_G_solar_a[i]; Q_G_solar_b[i]; Q_G_solar_c[i]] - [Q_Demand_a[i]; Q_Demand_b[i]; Q_Demand_c[i]] == imag([V_a_r[i] + im*(V_a_i[i]); V_b_r[i] + im*(V_b_i[i]); V_c_r[i] + im*(V_c_i[i])] .* conj([conj(Y[3*i-2,:]')*v_fase; conj(Y[3*i-1,:]')*v_fase; conj(Y[3*i,:]')*v_fase]) ))

# Model solution
JuMP.optimize!(model)  


########## 3. RESULTS OBTAINED ##########

# Extract the values of the solar generation variables
P_G_solar_a_val = value.(P_G_solar_a)
P_G_solar_b_val = value.(P_G_solar_b)
P_G_solar_c_val = value.(P_G_solar_c)
Q_G_solar_a_val = value.(Q_G_solar_a)
Q_G_solar_b_val = value.(Q_G_solar_b)
Q_G_solar_c_val = value.(Q_G_solar_c)

# Extract the values of the voltage variables (real and imaginary part)
V_a_r_val = value.(V_a_r)
V_a_i_val = value.(V_a_i)
V_b_r_val = value.(V_b_r)
V_b_i_val = value.(V_b_i)
V_c_r_val = value.(V_c_r)
V_c_i_val = value.(V_c_i)

# Calculate the magnitude and argument of the voltages (polar form)
V_mag_a = sqrt.(V_a_r_val.^2 .+ V_a_i_val.^2)  # Magnitude of phase A voltage
V_mag_b = sqrt.(V_b_r_val.^2 .+ V_b_i_val.^2)  # Magnitude of phase B voltage
V_mag_c = sqrt.(V_c_r_val.^2 .+ V_c_i_val.^2)  # Magnitude of phase C voltage
V_arg_a = rad2deg.(atan.(V_a_i_val, V_a_r_val))  # Argument of phase A voltage in degrees 
V_arg_b = rad2deg.(atan.(V_b_i_val, V_b_r_val))  # Argument of phase B voltage in degrees
V_arg_c = rad2deg.(atan.(V_c_i_val, V_c_r_val))  # Argument of phase C voltage in degrees

# Extract the values of the power variables
P_G_a_val = value.(P_G_a)
Q_G_a_val = value.(Q_G_a)
P_G_b_val = value.(P_G_b)
Q_G_b_val = value.(Q_G_b)
P_G_c_val = value.(P_G_c)
Q_G_c_val = value.(Q_G_c)

# Initialize the vectors to store the active power flow for each line
Pline_a_ij = zeros(nL)  # Empty array to store the active power flow in phase A
Pline_b_ij = zeros(nL)  # Empty array to store the active power flow in phase B
Pline_c_ij = zeros(nL)  # Empty array to store the active power flow in phase C

Pline_a_ji = zeros(nL)  # Empty array to store the active power flow in phase A
Pline_b_ji = zeros(nL)  # Empty array to store the active power flow in phase B
Pline_c_ji = zeros(nL)  # Empty array to store the active power flow in phase C

# Calculate the active power flow for each line
for k in 1:nL
    i = dLinea_fbus[k]
    j = dLinea_tbus[k]
    Yl=[(inv(Zij)+0.5.*Ysh)  -inv(Zij) ; -inv(Zij) (inv(Zij)+0.5.*Ysh)]
    v_val=[V_a_r_val[i] + im*(V_a_i_val[i]); V_b_r_val[i] + im*(V_b_i_val[i]); V_c_r_val[i] + im*(V_c_i_val[i]); V_a_r_val[j] + im*(V_a_i_val[j]); V_b_r_val[j] + im*(V_b_i_val[j]); V_c_r_val[j] + im*(V_c_i_val[j])]

    Pij=real([V_a_r_val[i] + im*(V_a_i_val[i]); V_b_r_val[i] + im*(V_b_i_val[i]); V_c_r_val[i] + im*(V_c_i_val[i])] .* conj([conj(Yl[1,:]')*v_val; conj(Yl[2,:]')*v_val; conj(Yl[3,:]')*v_val]))
    Pji=real([V_a_r_val[j] + im*(V_a_i_val[j]); V_b_r_val[j] + im*(V_b_i_val[j]); V_c_r_val[j] + im*(V_c_i_val[j])] .* conj([conj(Yl[4,:]')*v_val; conj(Yl[5,:]')*v_val; conj(Yl[6,:]')*v_val]))

    Pline_a_ij[k] = Pij[1]  # Store phase A power flow
    Pline_b_ij[k] = Pij[2]  # Store phase B power flow
    Pline_c_ij[k] = Pij[3]  # Store phase C power flow

    Pline_a_ji[k] = Pji[1]  # Store phase A power flow
    Pline_b_ji[k] = Pji[2]  # Store phase B power flow
    Pline_c_ji[k] = Pji[3]  # Store phase C power flow
end

# Initialize the vectors to store the reactive power flow for each line
Qline_a_ij = zeros(nL)  # Empty array to store the reactive power flow in phase A
Qline_b_ij = zeros(nL)  # Empty array to store the reactive power flow in phase B
Qline_c_ij = zeros(nL)  # Empty array to store the reactive power flow in phase C

Qline_a_ji = zeros(nL)  # Empty array to store the reactive power flow in phase A
Qline_b_ji = zeros(nL)  # Empty array to store the reactive power flow in phase B
Qline_c_ji = zeros(nL)  # Empty array to store the reactive power flow in phase C


# Calculate the reactive power flow for each line
for k in 1:nL
    i = dLinea_fbus[k]
    j = dLinea_tbus[k]
    Yl=[(inv(Zij)+0.5.*Ysh)  -inv(Zij) ; -inv(Zij) (inv(Zij)+0.5.*Ysh)]
    v_val=[V_a_r_val[i] + im*(V_a_i_val[i]); V_b_r_val[i] + im*(V_b_i_val[i]); V_c_r_val[i] + im*(V_c_i_val[i]); V_a_r_val[j] + im*(V_a_i_val[j]); V_b_r_val[j] + im*(V_b_i_val[j]); V_c_r_val[j] + im*(V_c_i_val[j])]

    Qij=imag([V_a_r_val[i] + im*(V_a_i_val[i]); V_b_r_val[i] + im*(V_b_i_val[i]); V_c_r_val[i] + im*(V_c_i_val[i])] .* conj([conj(Yl[1,:]')*v_val; conj(Yl[2,:]')*v_val; conj(Yl[3,:]')*v_val]))
    Qji=imag([V_a_r_val[j] + im*(V_a_i_val[j]); V_b_r_val[j] + im*(V_b_i_val[j]); V_c_r_val[j] + im*(V_c_i_val[j])] .* conj([conj(Yl[4,:]')*v_val; conj(Yl[5,:]')*v_val; conj(Yl[6,:]')*v_val]))

    Qline_a_ij[k] = Qij[1]  # Store phase A power flow
    Qline_b_ij[k] = Qij[2]  # Store phase B power flow
    Qline_c_ij[k] = Qij[3]  # Store phase C power flow

    Qline_a_ji[k] = Qji[1]  # Store phase A power flow
    Qline_b_ji[k] = Qji[2]  # Store phase B power flow
    Qline_c_ji[k] = Qji[3]  # Store phase C power flow


end

println("RESULTS")
println("Optimisation status: ", termination_status(model))
println("Value of the objective function:", objective_value(model))


# Loading the DataFrames package to organize and display the results in tables
using DataFrames

println("Nodal voltages in rectangular form (real and imaginary parts, in volts):")
df0 = DataFrame(
    Node = 1:nN,  
    V_a_r = V_a_r_val,  
    V_a_i = V_a_i_val, 
    V_b_r = V_b_r_val, 
    V_b_i = V_b_i_val, 
    V_c_r = V_c_r_val,  
    V_c_i = V_c_i_val   
)
println(df0)

println("Nodal voltages in polar form (magnitude in volts and angle in degrees):")
df1 = DataFrame(
    Node = 1:nN,  
    V_a_mag = V_mag_a,  
    V_a_arg = V_arg_a,  
    V_b_mag = V_mag_b,  
    V_b_arg = V_arg_b,  
    V_c_mag = V_mag_c,  
    V_c_arg = V_arg_c   
)
println(df1)

println("Active and reactive power generated by conventional generators (W / var):")
df2 = DataFrame(
    Nodo = 1:nN,  
    P_G_a = P_G_a_val,       
    Q_G_a = Q_G_a_val,       
    P_G_b = P_G_b_val,       
    Q_G_b = Q_G_b_val,       
    P_G_c = P_G_c_val,       
    Q_G_c = Q_G_c_val       
    )
println(df2)

println("Active and reactive power generated by photovoltaic generators (W / var):")
df3 = DataFrame(
    Node = 1:nN,        
    P_a_solar = P_G_solar_a_val,
    P_b_solar = P_G_solar_b_val,
    P_c_solar = P_G_solar_c_val,
    Q_a_solar = Q_G_solar_a_val,
    Q_b_solar = Q_G_solar_b_val,
    Q_c_solar = Q_G_solar_c_val
)
println(df3)

println("Active power flows in the lines (W):")
df4 = DataFrame(
    Line = 1:nL,        
    Initial_node = dLinea_fbus,  
    Final_node = dLinea_tbus,  
    Pij_a = Pline_a_ij,      
    Pij_b = Pline_b_ij,      
    Pij_c = Pline_c_ij,       
    Pji_a = Pline_a_ji,      
    Pji_b = Pline_b_ji,      
    Pji_c = Pline_c_ji     
)
println(df4)


println("Reactive power flows in the lines (var):")
df5 = DataFrame(
    Line = 1:nL,        
    Initial_node = dLinea_fbus,  
    Final_node = dLinea_tbus,  
    Qij_a = Qline_a_ij,      
    Qij_b = Qline_b_ij,      
    Qij_c = Qline_c_ij,      
    Qji_a = Qline_a_ji,      
    Qji_b = Qline_b_ji,      
    Qji_c = Qline_c_ji       
)
println(df5)