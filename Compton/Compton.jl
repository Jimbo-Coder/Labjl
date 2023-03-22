using Plots, DataFrames, LsqFit, LinearAlgebra

function dist()
   d = -[16.7, 15, 14, 13, 12, 11, 10, 8, 6, 4, 2] #centimeters
   Aexp = [2124., 1410., 1218., 1102., 1025., 987., 928., 898., 860., 840., 822.] #cps
   Aerr = [15.5, 13., 12., 11., 11, 11, 11, 11, 11, 11, 11]
   
   @.model(x,p) = p[1]* ((x - p[2])^(-2)) + p[3]

   fit1 = curve_fit(model, d, Aexp, Aerr, [3000.,-26.,800.])
   covar = estimate_covar(fit1)

   if (fit1.converged ==false)
     print("FAILED")  
   end

   g1=scatter(d,Aexp, yerror = 2*Aerr, label = "Experimental", dpi = 300,title = "Intensity vs Distance",xlabel = "Distance (cm)", ylabel = "Intensity (cps)")
   plot!(d, model(d,fit1.param), label = "Theor. Fit")
   gui()
   print("Fit Parameters ", fit1.param)
   print("\n")
   print("Errors ",sqrt.(LinearAlgebra.diag(covar)))
   
   g2=scatter(d,log.(Aexp), label = "log Experimental", dpi = 300,title = "log(Intensity) vs Distance",xlabel = "Distance (cm)", ylabel = "log(Intensity) (cps)")
   plot!(d, log.(model(d,fit1.param)), label = "log Theor. Fit")
   gui()

   #savefig(g1, "C:/Users/maxri/Desktop/Classes_4-2/phys382 Adv.Lab/Labjl/Compton/Comptondistance.png")
   #savefig(g2, "C:/Users/maxri/Desktop/Classes_4-2/phys382 Adv.Lab/Labjl/Compton/Comptondistancelog.png")
end


function lead()
   xl = [0., 0.6, 1.2, 1.8,2.3, 3.3, 4.05]
   Aexp = [1402, 662, 310.29, 145.07, 99.00, 46.0, 32.19]
   Aerr = [10., 10., 6., 4., 4., 2., 2.]

   @.model(x,p) = p[1]* (exp( - p[2]*x)) + p[3]

   fit1 = curve_fit(model, xl, Aexp, Aerr, [1500., 40., 20.])
   covar = estimate_covar(fit1)

   if (fit1.converged ==false)
     print("FAILED")  
   end

   g1=scatter(xl,Aexp, yerror = 2*Aerr, label = "Experimental", dpi = 300,title = "Intensity vs Lead Thickness",xlabel = "Lead Thickness (cm)", ylabel = "Intensity (cps)")
   plot!(xl, model(xl,fit1.param), label = "Theor. Fit")
   gui()
   print("Fit Parameters ", fit1.param)
   print("\n")
   print("Errors ",sqrt.(LinearAlgebra.diag(covar)))

   g2=scatter(xl,log.(Aexp), label = "log Experimental", dpi = 300,title = "log(Intensity) vs lead thickness",xlabel = "Distance (cm)", ylabel = "log(Intensity) (cps)")
   plot!(xl, log.(model(xl,fit1.param)), label = "log Theor. Fit")
   gui()

   #savefig(g1, "C:/Users/maxri/Desktop/Classes_4-2/phys382 Adv.Lab/Labjl/Compton/Comptonlead.png")
   #savefig(g2, "C:/Users/maxri/Desktop/Classes_4-2/phys382 Adv.Lab/Labjl/Compton/Comptondleadlog.png")
end


function compton()
   th = [0., 12.1, 25.9, 37., 48., 57.]
   Ep = [661.7, 623.34, 371.27, 660.86, 573.60,644.80]
   
   hc = (6.63e-34)*(3e8)
   Epj = Ep*(1000)*(1.60e-19)

   hmc = (6.63e-34)/ ((9.1e-31)*(3e8))

   λ = hc ./Epj
   λ0 = hc / Epj[1]
   dl = λ .- λ0

   @.model(x,p) = p[1]* (1 - cos(x*pi/180))

   fit1 = curve_fit(model, th, dl, [hmc])
   covar = estimate_covar(fit1)

   if (fit1.converged ==false)
      print("FAILED")  
   end
   
   g1=scatter(th,dl, label = "Experimental", dpi = 300,title = "Δλ vs Angle(θ)",xlabel = "Angle θ (deg)", ylabel = "Δλ (m)")
   plot!(th, model(th,fit1.param), label = "Theor. Fit")
   gui()
   print("Fit Parameters ", fit1.param)
   print("\n")
   print("Errors ",sqrt.(LinearAlgebra.diag(covar)))

   #savefig(g1, "C:/Users/maxri/Desktop/Classes_4-2/phys382 Adv.Lab/Labjl/Compton/Comptonscat.png")

end

