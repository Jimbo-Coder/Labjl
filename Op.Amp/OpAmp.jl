using Plots, DataFrames, LsqFit

function comparator()
   v1 = [5.,4.20,4.05,3.95,3.85,3.80,2.99,-5.00,-4.20,-4.06,-3.95,-3.85,-3.80,-3.00]
   v2 = 4*ones(7)
   append!(v2, -4*ones(7))
   v0 = [-7.63,-7.57,-7.50,8.51, 8.50, 8.50, 8.49, 8.50, 8.50, 8.50, -7.51, -7.56, -7.57, -7.63]

   p= Plots.scatter((v2-v1), v0, xlabel = "v2-v1(v)",ylabel = "V_out (v)", title = "comparator")
   gui()
   #Plots.savefig(p, "Comparatorplot.png")
end

function inverting()
   v1 = [1.970,1.532,.739,.514,.1035,-.1022,-.508,-.693,-1.502,-2.001]
   v2 = [-7.37,-7.29,-3.62,-2.51,-0.50,0.49,2.48,3.39,7.35,8.11]

   R2 = 1e4; R1 = 2e3; a = R2/R1;

   p= Plots.scatter(v1, v2, xlabel = "V_in(v)",ylabel = "V_out (v)", title = "Inverting Amplifier",label ="Experimental")
   Plots.plot!(v1, -a*v1, label = "Theoretical")
   #Plots.savefig(p, "InvertingAmp.png")
end


function instrumental()
   v1 = [4.172,4.226,4.260,4.265,4.276,4.282,4.290,4.290,4.327,4.375]
   v2 = 4.272*ones(length(v1))
   v0 = [8.50,8.50,6.14,5.42,3.90,3.00,0.9,0.43,-3.91,-7.57]
   
   R0 = 1e6; Rg = 2e4; A = (1 + 2*R0/Rg)
   vth = A *(v2-v1)

   p= Plots.scatter(v2-v1,v0, xlabel = "V2-V1(v)",ylabel = "V_out (v)", title = "Instrumental Amplifier",label ="Experimental")
   Plots.plot!(v2-v1, A*(v2-v1), label = "Theoretical")
   #Plots.savefig(p, "InstrumentalAmp.png")
end