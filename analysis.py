#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import ROOT as r
import numpy as np
import matplotlib.pyplot as plt
import math
from ROOT import RooFit
from ROOT import TChain
from ROOT import RooStats
from datetime import datetime
from uncertainties import ufloat
from scipy.optimize import curve_fit

plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

#get_ipython().run_line_magic('jsroot', 'on')

c = r.TCanvas()
r.gDirectory.cd(0)


# # Functions

# # Monte Carlo mass fits

# ## Monte Carlo mass fit, $B^+ \rightarrow \bar{D}^0\pi^+$

# In[ ]:


# import the data
tree = TChain("ST-b2oc")
tree.Add("more_more_data/combined/00266897_0000000*_1.highmult_2024-Friend-B2OC-W4042-DOWN.root") # path to MC data

if tree.GetEntries() == 0:
    print("TChain is empty. Please check the file paths.")
    exit()

# plotting the mass
# tree.Draw("Bp_DTF_OwnPV_MASS", "Bp_DTF_OwnPV_MASS < 5950", "logy")
# c.Draw()


# In[ ]:


# create a RooRealVar for the observable to fit
x = r.RooRealVar("Bp_DTF_OwnPV_MASS", "Bp_DTF_OwnPV_MASS", 4900, 5950)

# convert the data to RooDataSet
data = r.RooDataSet("data", "dataset from tree", r.RooArgSet(x), RooFit.Import(tree))

# create the pdfs
mean = r.RooRealVar("mean", "mean", 5278.46, 5250, 5300)

alpha = r.RooRealVar("alpha", "alpha", 1.25015, 1, 2)
n = r.RooRealVar("n", "n", 2.40748, 0.5, 5)
cb_sigma = r.RooRealVar("cb_sigma", "cb_sigma", 25.7691, 20, 30)
crystal_ball = r.RooCBShape("crystal_ball", "Crystal ball PDF", x, mean, cb_sigma, alpha, n)

tau = r.RooRealVar("tau", "Decay constant", -0.00219122, -1, 0)
background = r.RooExponential("background", "Exponential background", x, tau)

alpha_2 = r.RooRealVar("alpha2", "alpha2", -2.2701, -20, -0.01)
n_2 = r.RooRealVar("n2", "n2", 2.47972, 0.05, 50)
cb_sigma_2 = r.RooRealVar("cb_sigma2", "cb_sigma2", 12.7184, 5, 20)
crystal_ball_2 = r.RooCBShape("crystal_ball2", "Crystal ball PDF 2", x, mean, cb_sigma_2, alpha_2, n_2)

# combine the pdfs
frac_cb_2 = r.RooRealVar("frac_cb_2", "Fraction of crystal ball 2", 0.5, 0.0, 1.0)
frac_background = r.RooRealVar("frac_background", "Fraction of background", 0.2, 0.15, 0.25)
combined_pdf = r.RooAddPdf("combined_pdf", "Gaussian + crystal ball + background",
                           r.RooArgList(crystal_ball_2, crystal_ball, background),
                           r.RooArgList(frac_cb_2, frac_background))

# lock variables
mean.setConstant(True)
alpha.setConstant(False)
n.setConstant(False)
cb_sigma.setConstant(False)
tau.setConstant(True)
alpha_2.setConstant(False)
n_2.setConstant(False)
cb_sigma_2.setConstant(False)
frac_background.setConstant(False)

# fit to the data and plot
combined_pdf.fitTo(data, RooFit.PrintLevel(1))

frame = x.frame(RooFit.Title("2 crystal ball fit, B^{+} #rightarrow #bar{D}^{0}#pi^{+} Monte Carlo B^{+} mass"))
data.plotOn(frame)
combined_pdf.plotOn(frame)
combined_pdf.plotOn(frame, RooFit.Components(crystal_ball), RooFit.LineStyle(r.kDashed))
combined_pdf.plotOn(frame, RooFit.Components(crystal_ball_2), RooFit.LineStyle(r.kDashed), RooFit.LineColor(r.kRed))
combined_pdf.plotOn(frame, RooFit.Components(background), RooFit.LineStyle(r.kDashed), RooFit.LineColor(r.kGreen))

canvas = r.TCanvas("canvas", "2 crystal ball fit", 800, 600)
r.gDirectory.cd(0)
frame.Draw()
canvas.SetLogy()
frame.GetXaxis().SetTitle("Mass (MeV/c^{2})")
frame.SetMinimum(0.3)
canvas.Draw()


# In[ ]:


print(x, mean, cb_sigma, alpha, n, tau, alpha_2, n_2, cb_sigma_2, frac_cb_2, frac_background)


# ## Monte Carlo mass fit,  $B^+ \rightarrow J/\psi K^+$

# In[ ]:


# import the data
tree = TChain("ST-b2cc")
tree.Add("more_more_data/combined/00265775_0000000*_1.highmult_2024-Friend-B2CC-W4042-DOWN.root") # path to MC data

if tree.GetEntries() == 0:
    print("TChain is empty. Please check the file paths.")
    exit()

# plotting the mass
# tree.Draw("Bp_DTF_OwnPV_MASS", "Bp_DTF_OwnPV_MASS < 5750", "logy")
# c.Draw()


# In[ ]:


# create a RooRealVar for the observable to fit
x = r.RooRealVar("Bp_DTF_OwnPV_MASS", "Bp_DTF_OwnPV_MASS", 5050, 5550)

# convert the data to RooDataSet
data = r.RooDataSet("data", "dataset from tree", r.RooArgSet(x), RooFit.Import(tree))

# create the pdfs (2 crystal balls)
mean = r.RooRealVar("mean", "mean", 5279.38, 5250, 5300)

tau = r.RooRealVar("tau", "Decay constant", -0.00107319, -1, 0)
background = r.RooExponential("background", "Exponential background", x, tau)

alpha_1 = r.RooRealVar("alpha", "alpha", 1.59501, 0.1, 2)
n_1 = r.RooRealVar("n", "n", 2.75673, 0.5, 5)
cb_sigma_1 = r.RooRealVar("cb_sigma", "cb_sigma", 7.2753, 5, 10)
crystal_ball = r.RooCBShape("crystal_ball", "Crystal ball PDF", x, mean, cb_sigma_1, alpha_1, n_1)

alpha_2 = r.RooRealVar("alpha2", "alpha2", -1.29179, -1.5, -0.1)
n_2 = r.RooRealVar("n2", "n2", 5.29, 0.5, 10)
cb_sigma_2 = r.RooRealVar("cb_sigma2", "cb_sigma2", 8.21279, 5, 10)
crystal_ball_2 = r.RooCBShape("crystal_ball2", "Crystal ball PDF 2", x, mean, cb_sigma_2, alpha_2, n_2)

# combine the pdfs and fit the data
frac_cb_2 = r.RooRealVar("frac_cb_2", "Fraction of crystal ball 2", 0.5, 0.0, 1.0)
frac_background = r.RooRealVar("frac_background", "Fraction of background", 0.2, 0.0, 1.0)
combined_pdf = r.RooAddPdf("combined_pdf", "Gaussian + crystal ball + background",
                           r.RooArgList(crystal_ball_2, crystal_ball, background),
                           r.RooArgList(frac_cb_2, frac_background))

mean.setConstant(True)
tau.setConstant(True)
alpha_1.setConstant(True)
n_1.setConstant(True)
cb_sigma_1.setConstant(True)
alpha_2.setConstant(True)
n_2.setConstant(True)
cb_sigma_2.setConstant(True)

combined_pdf.fitTo(data, RooFit.PrintLevel(1))

# plot it
frame = x.frame(RooFit.Title("2 crystal ball fit, B^{+} #rightarrow J/#psiK^{+} Monte Carlo B^{+} mass"))
data.plotOn(frame)
combined_pdf.plotOn(frame)
combined_pdf.plotOn(frame, RooFit.Components(crystal_ball), RooFit.LineStyle(r.kDashed))
combined_pdf.plotOn(frame, RooFit.Components(crystal_ball_2), RooFit.LineStyle(r.kDashed), RooFit.LineColor(r.kRed))
combined_pdf.plotOn(frame, RooFit.Components(background), RooFit.LineStyle(r.kDashed), RooFit.LineColor(r.kGreen))

canvas = r.TCanvas("canvas", "2x crystal ball fit", 800, 600)
r.gDirectory.cd(0)
frame.Draw()
canvas.SetLogy()
frame.SetMinimum(100)
frame.GetXaxis().SetTitle("Mass (MeV/c^{2})")
canvas.Draw()


# In[ ]:


print(x, mean, cb_sigma_1, alpha_1, n_1, cb_sigma_2, alpha_2, n_2, tau, frac_cb_2, frac_background)


# # Time binned mass, $B^+ \rightarrow \bar{D}^0\pi^+$

# ## Up

# In[ ]:


# import the data
tree_mdpi_up = TChain("ST-b2oc")
tree_mdpi_up.Add("more_more_data/combined/0028*_1.highstats-Small-B2OC-UP.root") # path to data

if tree_mdpi_up.GetEntries() == 0:
    print("TChain is empty. Please check the file paths.")
    exit()


# In[ ]:


# create a RooRealVar for the observable to fit
x = r.RooRealVar("Bp_M", "Bp_M", 5200, 5950)

# this to have time in the dataset for later
min_time_dpi_up = tree_mdpi_up.GetMinimum("GPSTIME")
max_time_dpi_up = tree_mdpi_up.GetMaximum("GPSTIME")
print(min_time_dpi_up, max_time_dpi_up)
print(datetime.fromtimestamp(min_time_dpi_up / 1e6), datetime.fromtimestamp(max_time_dpi_up / 1e6))
time = r.RooRealVar("GPSTIME", "GPSTIME", min_time_dpi_up, max_time_dpi_up)
block = r.RooRealVar("block", "block", 1, 8)

# convert the data to RooDataSet
data_dpi_up = r.RooDataSet("data_dpi_up", "dataset from tree", r.RooArgSet(x, time, block), RooFit.Import(tree_mdpi_up))


# In[ ]:


data_dpi_up_b5 = data_dpi_up.reduce("block >= 5 && block < 6")
data_dpi_up_b8 = data_dpi_up.reduce("block >= 8 && block < 9")


# In[ ]:


print(data_dpi_up_b5.numEntries())
print(data_dpi_up_b8.numEntries())


# ### Plots

# In[ ]:


c = r.TCanvas()
r.gDirectory.cd(0)

tree_mdpi_up.Draw("block >> histblock(8, 1, 9)")

c.Draw()


# In[ ]:


c = r.TCanvas()
r.gDirectory.cd(0)

tree_mdpi_up.Draw("Bp_M >> histmass(100, 4875, 5900)")

hist = r.gPad.GetPrimitive("histmass")

hist.GetXaxis().SetTitle("Mass (MeV/c^{2})")
hist.GetYaxis().SetTitle("Events")
hist.SetTitle(f"B^{{+}} mass, B^{{+}} #rightarrow #bar{{D}}^{{0}}#pi^{{+}}")

hist.SetStats(False)

r.gPad.Update()

c.Draw()


# ### Block 5

# In[ ]:


# extract and sort all the times
n_entries = data_dpi_up_b5.numEntries()
times_dpi_up_b5 = np.empty(n_entries, dtype=np.float64)

for i in range(n_entries):
    values = data_dpi_up_b5.get(i)
    times_dpi_up_b5[i] = values.getRealValue("GPSTIME")

times_dpi_up_b5.sort()

# define bin edges
n_bins = 10
time_bins_up_b5 = np.quantile(times_dpi_up_b5, np.linspace(0, 1, n_bins+1))

print(time_bins_up_b5)


# In[ ]:


signal_yields_dpi_up_b5 = np.empty(10, dtype=object)

background_yields_dpi_up_b5 = np.empty(10, dtype=object)


# In[ ]:


for i in range(n_bins):
    time_min = time_bins_up_b5[i]
    time_max = time_bins_up_b5[i+1]

    # define the cut string for this time bin
    cut_str = f"GPSTIME >= {time_min} && GPSTIME < {time_max}"

    # create the binned dataset
    binned_data_dpi_up = data_dpi_up_b5.reduce(cut_str)

    print(f"Bin {i+1}: {time_min} - {time_max}, entries: {binned_data_dpi_up.numEntries()}")

    # create the pdfs
    mean = r.RooRealVar("mean", "mean", 5277.8, 5277, 5280)

    cb_sigma_oldg = r.RooRealVar("cb_sigma_oldg", "cb_sigma_oldg", 18.39, 15, 25)
    alpha_oldg = r.RooRealVar("alpha_oldg", "alpha_oldg", -1.42686, -2, -0.8)
    n_oldg = r.RooRealVar("n_oldg", "n_oldg", 19.9954, 10, 20)
    cb_oldg = r.RooCBShape("cb_oldg", "Second CB PDF", x, mean, cb_sigma_oldg, alpha_oldg, n_oldg)

    alpha = r.RooRealVar("alpha", "alpha", 0.6, 0.5, 1)
    n = r.RooRealVar("n", "n", 19.6699, 15, 25)
    cb_sigma = r.RooRealVar("cb_sigma", "cb_sigma", 24.7, 20, 50)
    crystal_ball = r.RooCBShape("crystal_ball", "Crystal ball PDF", x, mean, cb_sigma, alpha, n)

    tau = r.RooRealVar("tau", "Decay constant", -0.00051, -0.0006, 0)
    background = r.RooExponential("background", "Exponential background", x, tau)

    mean.setConstant(True)
    tau.setConstant(True)
    cb_sigma_oldg.setConstant(True)
    alpha_oldg.setConstant(True)
    n_oldg.setConstant(True)
    alpha.setConstant(True) #
    n.setConstant(True)
    cb_sigma.setConstant(True)

    # use yields instead of fractions
    ntot = r.RooRealVar("ntot", "total yield", binned_data_dpi_up.numEntries(), 0, data_dpi_up.numEntries() * 10)

#    ncb1 = r.RooRealVar("ncb1", "yield cb1", binned_data_dpi_up.numEntries() / 10, binned_data_dpi_up.numEntries() / 40, binned_data_dpi_up.numEntries())
#    ncb2 = r.RooRealVar("ncb2", "yield cb2", binned_data_dpi_up.numEntries() / 10, binned_data_dpi_up.numEntries() / 40, binned_data_dpi_up.numEntries())
    nsig = r.RooRealVar("nsig", "yield sig", binned_data_dpi_up.numEntries() * 0.15, binned_data_dpi_up.numEntries() * 0.12, binned_data_dpi_up.numEntries() * 0.4)
    nbkg = r.RooRealVar("nbkg", "yield bkg", binned_data_dpi_up.numEntries() * 0.85, binned_data_dpi_up.numEntries() * 0.6, binned_data_dpi_up.numEntries() * 0.88)
    f_cb1 = r.RooRealVar("f_cb1", "fraction of cb1", 0.6)

    # build the model and fit the data

    cb_sum = r.RooAddPdf("cb_sum", "signal pdf", r.RooArgList(cb_oldg, crystal_ball), r.RooArgList(f_cb1))

    combined_pdf_dpi_up_timebinned = r.RooAddPdf("combined_pdf", "2x crystal ball + background",
                               r.RooArgList(cb_sum, background),
                               r.RooArgList(nsig, nbkg))

    # fit the model
    fit_result = combined_pdf_dpi_up_timebinned.fitTo(binned_data_dpi_up,
                                                      #RooFit.Minos(True), # asymmetric uncertainties
                                                      #RooFit.Extended(True),
                                                      RooFit.Save(),
                                                      Strategy=2)

    print(x, mean, cb_sigma_oldg, alpha_oldg, n_oldg, cb_sigma, alpha, n, tau, ntot, nsig, nbkg)
    print(f"Background ratio: {nbkg.getVal() / ntot.getVal()}")
#    print(f"CB1/nsig ratio: {ncb1.getVal() / (ncb1.getVal() + ncb2.getVal())}")
    print(f"Status: {fit_result.status()}, CovQual: {fit_result.covQual()}")
    fit_result.Print()

    # extract the yield
    #nsig = r.RooAddition("nsig", "ncb1 + ncb2", r.RooArgList(ncb1, ncb2))

    #ncb1_err = ufloat(ncb1.getVal(), ncb1.getError())
    #ncb2_err = ufloat(ncb2.getVal(), ncb2.getError())

    #signal_yield = ncb1_err + ncb2_err

    signal_yield = ufloat(nsig.getVal(), nsig.getError())
    background_yield = ufloat(nbkg.getVal(), nbkg.getError())

    signal_yields_dpi_up_b5[i] = signal_yield
    background_yields_dpi_up_b5[i] = background_yield

    #frame = x.frame(RooFit.Title(f"Block 5 up Dpi, time bin {i+1}: {datetime.fromtimestamp(time_min / 1e6)} - {datetime.fromtimestamp(time_max / 1e6)}"))
    frame = x.frame(RooFit.Title(f"B^{{+}} mass fit, B^{{+}} #rightarrow #bar{{D}}^{{0}}#pi^{{+}}, block 5 (up), time bin {i+1}: {datetime.fromtimestamp(time_min / 1e6).date()} - {datetime.fromtimestamp(time_max / 1e6).date()}"))
    binned_data_dpi_up.plotOn(frame)
    combined_pdf_dpi_up_timebinned.plotOn(frame)
    combined_pdf_dpi_up_timebinned.plotOn(frame, RooFit.Components(crystal_ball), RooFit.LineStyle(r.kDashed))
    combined_pdf_dpi_up_timebinned.plotOn(frame, RooFit.Components(cb_oldg), RooFit.LineStyle(r.kDashed), RooFit.LineColor(r.kRed))
    combined_pdf_dpi_up_timebinned.plotOn(frame, RooFit.Components(background), RooFit.LineStyle(r.kDashed), RooFit.LineColor(r.kGreen))

    canvas = r.TCanvas(f"c_bin_{i+1}", f"Fit for time bin {i+1}", 800, 600)
    frame.GetXaxis().SetTitle("Mass (MeV/c^{2})")
    frame.Draw()
    frame.SetMinimum(0)
    frame.GetXaxis().SetRangeUser(5200, 5950)
    canvas.SaveAs(f"fulldata_b5_up_time_binned_dpi_fit_bin_{i+1}.png")
    frame.Draw()
    canvas.SetLogy()
    frame.SetMinimum(4000)
    canvas.SaveAs(f"fulldata_b5_up_time_binned_dpi_fit_bin_{i+1}_log.png")

    print(f"Yield for {datetime.fromtimestamp(time_min / 1e6)} - {datetime.fromtimestamp(time_max / 1e6)}: {signal_yield:.2f}\n\n")# ± {signal_error:.2f}")

print(f"All yields: {signal_yields_dpi_up_b5}")
print(f"Background yields: {background_yields_dpi_up_b5}")


# In[ ]:


# compute sWeights if needed
# splot_dpi_up = RooStats.SPlot("splot_dpi_b5_up", "sPlot Dpi block 5 up", data_dpi_up_b5, combined_pdf_dpi_up_timebinned, r.RooArgList(nsig, nbkg))


# In[ ]:


# print("sWeight branches:", [var.GetName() for var in data_dpi_up_b5.get() if var.GetName().endswith("_sw")])


# ### Block 8

# In[ ]:


# extract and sort all the times
n_entries = data_dpi_up_b8.numEntries()
times_dpi_up_b8 = np.empty(n_entries, dtype=np.float64)

for i in range(n_entries):
    values = data_dpi_up_b8.get(i)
    times_dpi_up_b8[i] = values.getRealValue("GPSTIME")

times_dpi_up_b8.sort()

# define bin edges
n_bins = 10
time_bins_up_b8 = np.quantile(times_dpi_up_b8, np.linspace(0, 1, n_bins+1))

print(time_bins_up_b8)


# In[ ]:


signal_yields_dpi_up_b8 = np.empty(10, dtype=object)

background_yields_dpi_up_b8 = np.empty(10, dtype=object)


# In[ ]:


for i in range(n_bins):
    time_min = time_bins_up_b8[i]
    time_max = time_bins_up_b8[i+1]

    # define the cut string for this time bin
    cut_str = f"GPSTIME >= {time_min} && GPSTIME < {time_max}"

    # create the binned dataset
    binned_data_dpi_up = data_dpi_up_b8.reduce(cut_str)

    print(f"Bin {i+1}: {time_min} - {time_max}, entries: {binned_data_dpi_up.numEntries()}")

    # create the pdfs
    mean = r.RooRealVar("mean", "mean", 5277.8, 5277, 5280)

    cb_sigma_oldg = r.RooRealVar("cb_sigma_oldg", "cb_sigma_oldg", 18.39, 15, 25)
    alpha_oldg = r.RooRealVar("alpha_oldg", "alpha_oldg", -1.42686, -2, -0.8)
    n_oldg = r.RooRealVar("n_oldg", "n_oldg", 19.9954, 10, 20)
    cb_oldg = r.RooCBShape("cb_oldg", "Second CB PDF", x, mean, cb_sigma_oldg, alpha_oldg, n_oldg)

    alpha = r.RooRealVar("alpha", "alpha", 0.618402, 0.4, 1)
    n = r.RooRealVar("n", "n", 5, 4, 15)
    cb_sigma = r.RooRealVar("cb_sigma", "cb_sigma", 24.7, 20, 50)
    crystal_ball = r.RooCBShape("crystal_ball", "Crystal ball PDF", x, mean, cb_sigma, alpha, n)

    tau = r.RooRealVar("tau", "Decay constant", -0.000395081, -0.0006, 0)
    background = r.RooExponential("background", "Exponential background", x, tau)

    mean.setConstant(True)
    tau.setConstant(True)
    cb_sigma_oldg.setConstant(True)
    alpha_oldg.setConstant(True)
    n_oldg.setConstant(True)
    alpha.setConstant(True) #
    n.setConstant(True)
    cb_sigma.setConstant(True)

    # use yields instead of fractions
    ntot = r.RooRealVar("ntot", "total yield", binned_data_dpi_up.numEntries(), 0, data_dpi_up.numEntries() * 10)

#    ncb1 = r.RooRealVar("ncb1", "yield cb1", binned_data_dpi_up.numEntries() / 10, binned_data_dpi_up.numEntries() / 40, binned_data_dpi_up.numEntries())
#    ncb2 = r.RooRealVar("ncb2", "yield cb2", binned_data_dpi_up.numEntries() / 10, binned_data_dpi_up.numEntries() / 40, binned_data_dpi_up.numEntries())
    nsig = r.RooRealVar("nsig", "yield sig", binned_data_dpi_up.numEntries() * 0.15, binned_data_dpi_up.numEntries() * 0.1, binned_data_dpi_up.numEntries() * 0.17)
    nbkg = r.RooRealVar("nbkg", "yield bkg", binned_data_dpi_up.numEntries() * 0.85, binned_data_dpi_up.numEntries() * 0.81, binned_data_dpi_up.numEntries() * 0.91)
    f_cb1 = r.RooRealVar("f_cb1", "fraction of cb1", 0.6)

    # build the model and fit the data

    cb_sum = r.RooAddPdf("cb_sum", "signal pdf", r.RooArgList(cb_oldg, crystal_ball), r.RooArgList(f_cb1))

    combined_pdf_dpi_up_timebinned = r.RooAddPdf("combined_pdf", "2x crystal ball + background",
                               r.RooArgList(cb_sum, background),
                               r.RooArgList(nsig, nbkg))

    # fit the model
    fit_result = combined_pdf_dpi_up_timebinned.fitTo(binned_data_dpi_up,
                                                      #RooFit.Minos(True), # asymmetric uncertainties
                                                      #RooFit.Extended(True),
                                                      RooFit.Save(),
                                                      Strategy=2)

    print(x, mean, cb_sigma_oldg, alpha_oldg, n_oldg, cb_sigma, alpha, n, tau, ntot, nsig, nbkg)
    print(f"Background ratio: {nbkg.getVal() / ntot.getVal()}")
#    print(f"CB1/nsig ratio: {ncb1.getVal() / (ncb1.getVal() + ncb2.getVal())}")
    print(f"Status: {fit_result.status()}, CovQual: {fit_result.covQual()}")
    fit_result.Print()

    # extract the yield
    #nsig = r.RooAddition("nsig", "ncb1 + ncb2", r.RooArgList(ncb1, ncb2))

    #ncb1_err = ufloat(ncb1.getVal(), ncb1.getError())
    #ncb2_err = ufloat(ncb2.getVal(), ncb2.getError())

    #signal_yield = ncb1_err + ncb2_err

    signal_yield = ufloat(nsig.getVal(), nsig.getError())
    background_yield = ufloat(nbkg.getVal(), nbkg.getError())

    signal_yields_dpi_up_b8[i] = signal_yield
    background_yields_dpi_up_b8[i] = background_yield

    frame = x.frame(RooFit.Title(f"B^{{+}} mass fit, B^{{+}} #rightarrow #bar{{D}}^{{0}}#pi^{{+}}, block 8 (up), time bin {i+1}: {datetime.fromtimestamp(time_min / 1e6).date()} - {datetime.fromtimestamp(time_max / 1e6).date()}"))
    binned_data_dpi_up.plotOn(frame)
    combined_pdf_dpi_up_timebinned.plotOn(frame)
    combined_pdf_dpi_up_timebinned.plotOn(frame, RooFit.Components(crystal_ball), RooFit.LineStyle(r.kDashed))
    combined_pdf_dpi_up_timebinned.plotOn(frame, RooFit.Components(cb_oldg), RooFit.LineStyle(r.kDashed), RooFit.LineColor(r.kRed))
    combined_pdf_dpi_up_timebinned.plotOn(frame, RooFit.Components(background), RooFit.LineStyle(r.kDashed), RooFit.LineColor(r.kGreen))

    canvas = r.TCanvas(f"c_bin_{i+1}", f"Fit for time bin {i+1}", 800, 600)
    frame.GetXaxis().SetTitle("Mass (MeV/c^{2})")
    frame.Draw()
    frame.SetMinimum(0)
    canvas.SaveAs(f"fulldata_b8_up_time_binned_dpi_fit_bin_{i+1}.png")
    frame.Draw()
    canvas.SetLogy()
    frame.SetMaximum(10001)
    frame.SetMinimum(2000)
    canvas.SaveAs(f"fulldata_b8_up_time_binned_dpi_fit_bin_{i+1}_log.png")

    print(f"Yield for {datetime.fromtimestamp(time_min / 1e6)} - {datetime.fromtimestamp(time_max / 1e6)}: {signal_yield:.2f}\n\n")# ± {signal_error:.2f}")

print(f"All yields: {signal_yields_dpi_up_b8}")
print(f"Background yields: {background_yields_dpi_up_b8}")


# ## Down

# In[ ]:


# import the data
tree_mdpi_down = TChain("ST-b2oc")
tree_mdpi_down.Add("more_more_data/combined/0028*_1.highstats-Small-B2OC-DOWN.root") # path to data

if tree_mdpi_down.GetEntries() == 0:
    print("TChain is empty. Please check the file paths.")
    exit()


# In[ ]:


# create a RooRealVar for the observable to fit
x = r.RooRealVar("Bp_M", "Bp_M", 5200, 5950) # 5950

# this to have time in the dataset for later
min_time_dpi_down = tree_mdpi_down.GetMinimum("GPSTIME")
max_time_dpi_down = tree_mdpi_down.GetMaximum("GPSTIME")
print(min_time_dpi_down, max_time_dpi_down)
print(datetime.fromtimestamp(min_time_dpi_down / 1e6), datetime.fromtimestamp(max_time_dpi_down / 1e6))
time = r.RooRealVar("GPSTIME", "GPSTIME", min_time_dpi_down, max_time_dpi_down)
block = r.RooRealVar("block", "block", 1, 8)

# convert the data to RooDataSet
data_dpi_down = r.RooDataSet("data_dpi_up", "dataset from tree", r.RooArgSet(x, time, block), RooFit.Import(tree_mdpi_down))


# In[ ]:


data_dpi_down_b6 = data_dpi_down.reduce("block >= 6 && block < 7")
data_dpi_down_b7 = data_dpi_down.reduce("block >= 7 && block < 8")


# In[ ]:


print(data_dpi_down_b6.numEntries())
print(data_dpi_down_b7.numEntries())


# ### Block 6

# In[ ]:


# extract and sort all the times
n_entries = data_dpi_down_b6.numEntries()
times_dpi_down_b6 = np.empty(n_entries, dtype=np.float64)

for i in range(n_entries):
    values = data_dpi_down_b6.get(i)
    times_dpi_down_b6[i] = values.getRealValue("GPSTIME")

times_dpi_down_b6.sort()

# define bin edges
n_bins = 10
time_bins_down_b6 = np.quantile(times_dpi_down_b6, np.linspace(0, 1, n_bins+1))

print(time_bins_down_b6)


# In[ ]:


signal_yields_dpi_down_b6 = np.empty(10, dtype=object)

background_yields_dpi_down_b6 = np.empty(10, dtype=object)


# In[ ]:


for i in range(n_bins):
    time_min = time_bins_down_b6[i]
    time_max = time_bins_down_b6[i+1]

    # define the cut string for this time bin
    cut_str = f"GPSTIME >= {time_min} && GPSTIME < {time_max}"

    # create the binned dataset
    binned_data_dpi_down = data_dpi_down_b6.reduce(cut_str)

    print(f"Bin {i+1}: {time_min} - {time_max}, entries: {binned_data_dpi_down.numEntries()}")

    # create the pdfs
    mean = r.RooRealVar("mean", "mean", 5277.8, 5277, 5280)

    cb_sigma_oldg = r.RooRealVar("cb_sigma_oldg", "cb_sigma_oldg", 18.39, 15, 25)
    alpha_oldg = r.RooRealVar("alpha_oldg", "alpha_oldg", -1.42686, -2, -0.8)
    n_oldg = r.RooRealVar("n_oldg", "n_oldg", 19.9954, 10, 20)
    cb_oldg = r.RooCBShape("cb_oldg", "Second CB PDF", x, mean, cb_sigma_oldg, alpha_oldg, n_oldg)

    alpha = r.RooRealVar("alpha", "alpha", 0.6, 0.4, 1)
    n = r.RooRealVar("n", "n", 19.6699, 15, 25)
    cb_sigma = r.RooRealVar("cb_sigma", "cb_sigma", 24.7, 20, 50)
    crystal_ball = r.RooCBShape("crystal_ball", "Crystal ball PDF", x, mean, cb_sigma, alpha, n)

    tau = r.RooRealVar("tau", "Decay constant", -0.00051, -0.0006, 0)
    background = r.RooExponential("background", "Exponential background", x, tau)

    mean.setConstant(True)
    tau.setConstant(True)
    cb_sigma_oldg.setConstant(True)
    alpha_oldg.setConstant(True)
    n_oldg.setConstant(True)
    alpha.setConstant(True) #
    n.setConstant(True)
    cb_sigma.setConstant(True)

    # use yields instead of fractions
    ntot = r.RooRealVar("ntot", "total yield", binned_data_dpi_down.numEntries(), 0, data_dpi_down.numEntries() * 10)

    nsig = r.RooRealVar("nsig", "yield sig", binned_data_dpi_down.numEntries() * 0.15, binned_data_dpi_down.numEntries() * 0.12, binned_data_dpi_down.numEntries() * 0.2)
    nbkg = r.RooRealVar("nbkg", "yield bkg", binned_data_dpi_down.numEntries() * 0.85, binned_data_dpi_down.numEntries() * 0.81, binned_data_dpi_down.numEntries() * 0.88)
    f_cb1 = r.RooRealVar("f_cb1", "fraction of cb1", 0.6)

    # build the model and fit the data

    cb_sum = r.RooAddPdf("cb_sum", "signal pdf", r.RooArgList(cb_oldg, crystal_ball), r.RooArgList(f_cb1))

    combined_pdf_dpi_down_timebinned = r.RooAddPdf("combined_pdf", "2x crystal ball + background",
                               r.RooArgList(cb_sum, background),
                               r.RooArgList(nsig, nbkg))

    # fit the model
    fit_result = combined_pdf_dpi_down_timebinned.fitTo(binned_data_dpi_down,
                                                      #RooFit.Minos(True), # asymmetric uncertainties
                                                      #RooFit.Extended(True),
                                                      RooFit.Save(),
                                                      Strategy=2)

    print(x, mean, cb_sigma_oldg, alpha_oldg, n_oldg, cb_sigma, alpha, n, tau, ntot, nsig, nbkg)
    print(f"Background ratio: {nbkg.getVal() / ntot.getVal()}")
#    print(f"CB1/nsig ratio: {ncb1.getVal() / (ncb1.getVal() + ncb2.getVal())}")
    print(f"Status: {fit_result.status()}, CovQual: {fit_result.covQual()}")
    fit_result.Print()

    signal_yield = ufloat(nsig.getVal(), nsig.getError())
    background_yield = ufloat(nbkg.getVal(), nbkg.getError())

    signal_yields_dpi_down_b6[i] = signal_yield
    background_yields_dpi_down_b6[i] = background_yield

    frame = x.frame(RooFit.Title(f"B^{{+}} mass fit, B^{{+}} #rightarrow #bar{{D}}^{{0}}#pi^{{+}}, block 6 (down), time bin {i+1}: {datetime.fromtimestamp(time_min / 1e6).date()} - {datetime.fromtimestamp(time_max / 1e6).date()}"))
    binned_data_dpi_down.plotOn(frame)
    combined_pdf_dpi_down_timebinned.plotOn(frame)
    combined_pdf_dpi_down_timebinned.plotOn(frame, RooFit.Components(crystal_ball), RooFit.LineStyle(r.kDashed))
    combined_pdf_dpi_down_timebinned.plotOn(frame, RooFit.Components(cb_oldg), RooFit.LineStyle(r.kDashed), RooFit.LineColor(r.kRed))
    combined_pdf_dpi_down_timebinned.plotOn(frame, RooFit.Components(background), RooFit.LineStyle(r.kDashed), RooFit.LineColor(r.kGreen))

    canvas = r.TCanvas(f"c_bin_{i+1}", f"Fit for time bin {i+1}", 800, 600)
    frame.GetXaxis().SetTitle("Mass (MeV/c^{2})")
    frame.Draw()
    frame.SetMinimum(0)
    canvas.SaveAs(f"fulldata_b6_down_time_binned_dpi_fit_bin_{i+1}.png")
    frame.Draw()
    canvas.SetLogy()
    frame.SetMinimum(3000)
    canvas.SaveAs(f"fulldata_b6_down_time_binned_dpi_fit_bin_{i+1}_log.png")

    print(f"Yield for {datetime.fromtimestamp(time_min / 1e6)} - {datetime.fromtimestamp(time_max / 1e6)}: {signal_yield:.2f}\n\n")# ± {signal_error:.2f}")

print(f"All yields: {signal_yields_dpi_down_b6}")
print(f"Background yields: {background_yields_dpi_down_b6}")


# ### Block 7

# In[ ]:


# extract and sort all the times
n_entries = data_dpi_down_b7.numEntries()
times_dpi_down_b7 = np.empty(n_entries, dtype=np.float64)

for i in range(n_entries):
    values = data_dpi_down_b7.get(i)
    times_dpi_down_b7[i] = values.getRealValue("GPSTIME")

times_dpi_down_b7.sort()

# define bin edges
n_bins = 10
time_bins_down_b7 = np.quantile(times_dpi_down_b7, np.linspace(0, 1, n_bins+1))

print(time_bins_down_b7)


# In[ ]:


signal_yields_dpi_down_b7 = np.empty(10, dtype=object)

background_yields_dpi_down_b7 = np.empty(10, dtype=object)


# In[ ]:


for i in range(n_bins):
    time_min = time_bins_down_b7[i]
    time_max = time_bins_down_b7[i+1]

    # define the cut string for this time bin
    cut_str = f"GPSTIME >= {time_min} && GPSTIME < {time_max}"

    # create the binned dataset
    binned_data_dpi_down = data_dpi_down_b7.reduce(cut_str)

    print(f"Bin {i+1}: {time_min} - {time_max}, entries: {binned_data_dpi_down.numEntries()}")

    # create the pdfs
    mean = r.RooRealVar("mean", "mean", 5277.8, 5277, 5280)

    cb_sigma_oldg = r.RooRealVar("cb_sigma_oldg", "cb_sigma_oldg", 18.39, 15, 25)
    alpha_oldg = r.RooRealVar("alpha_oldg", "alpha_oldg", -1.42686, -2, -0.8)
    n_oldg = r.RooRealVar("n_oldg", "n_oldg", 19.9954, 10, 20)
    cb_oldg = r.RooCBShape("cb_oldg", "Second CB PDF", x, mean, cb_sigma_oldg, alpha_oldg, n_oldg)

    alpha = r.RooRealVar("alpha", "alpha", 0.6, 0.3, 1)
    n = r.RooRealVar("n", "n", 7, 0.01, 10)
    cb_sigma = r.RooRealVar("cb_sigma", "cb_sigma", 24, 20, 30)
    crystal_ball = r.RooCBShape("crystal_ball", "Crystal ball PDF", x, mean, cb_sigma, alpha, n)

    tau = r.RooRealVar("tau", "Decay constant", -0.000411144, -0.0006, 0)
    background = r.RooExponential("background", "Exponential background", x, tau)

    mean.setConstant(True)
    tau.setConstant(True)
    cb_sigma_oldg.setConstant(True)
    alpha_oldg.setConstant(True)
    n_oldg.setConstant(True)
    alpha.setConstant(True) #
    n.setConstant(True)
    cb_sigma.setConstant(True)

    # use yields instead of fractions
    ntot = r.RooRealVar("ntot", "total yield", binned_data_dpi_down.numEntries(), 0, data_dpi_down.numEntries() * 10)

    nsig = r.RooRealVar("nsig", "yield sig", binned_data_dpi_down.numEntries() * 0.15, binned_data_dpi_down.numEntries() * 0.12, binned_data_dpi_down.numEntries() * 0.17)
    nbkg = r.RooRealVar("nbkg", "yield bkg", binned_data_dpi_down.numEntries() * 0.85, binned_data_dpi_down.numEntries() * 0.81, binned_data_dpi_down.numEntries() * 0.88)
    f_cb1 = r.RooRealVar("f_cb1", "fraction of cb1", 0.6)

    # build the model and fit the data

    cb_sum = r.RooAddPdf("cb_sum", "signal pdf", r.RooArgList(cb_oldg, crystal_ball), r.RooArgList(f_cb1))

    combined_pdf_dpi_down_timebinned = r.RooAddPdf("combined_pdf", "2x crystal ball + background",
                               r.RooArgList(cb_sum, background),
                               r.RooArgList(nsig, nbkg))

    # fit the model
    fit_result = combined_pdf_dpi_down_timebinned.fitTo(binned_data_dpi_down,
                                                      #RooFit.Minos(True), # asymmetric uncertainties
                                                      #RooFit.Extended(True),
                                                      RooFit.Save(),
                                                      Strategy=2)

    print(x, mean, cb_sigma_oldg, alpha_oldg, n_oldg, cb_sigma, alpha, n, tau, ntot, nsig, nbkg)
    print(f"Background ratio: {nbkg.getVal() / ntot.getVal()}")
#    print(f"CB1/nsig ratio: {ncb1.getVal() / (ncb1.getVal() + ncb2.getVal())}")
    print(f"Status: {fit_result.status()}, CovQual: {fit_result.covQual()}")
    fit_result.Print()

    signal_yield = ufloat(nsig.getVal(), nsig.getError())
    background_yield = ufloat(nbkg.getVal(), nbkg.getError())

    signal_yields_dpi_down_b7[i] = signal_yield
    background_yields_dpi_down_b7[i] = background_yield

    frame = x.frame(RooFit.Title(f"B^{{+}} mass fit, B^{{+}} #rightarrow #bar{{D}}^{{0}}#pi^{{+}}, block 7 (down), time bin {i+1}: {datetime.fromtimestamp(time_min / 1e6).date()} - {datetime.fromtimestamp(time_max / 1e6).date()}"))
    binned_data_dpi_down.plotOn(frame)
    combined_pdf_dpi_down_timebinned.plotOn(frame)
    combined_pdf_dpi_down_timebinned.plotOn(frame, RooFit.Components(crystal_ball), RooFit.LineStyle(r.kDashed))
    combined_pdf_dpi_down_timebinned.plotOn(frame, RooFit.Components(cb_oldg), RooFit.LineStyle(r.kDashed), RooFit.LineColor(r.kRed))
    combined_pdf_dpi_down_timebinned.plotOn(frame, RooFit.Components(background), RooFit.LineStyle(r.kDashed), RooFit.LineColor(r.kGreen))

    canvas = r.TCanvas(f"c_bin_{i+1}", f"Fit for time bin {i+1}", 800, 600)
    frame.GetXaxis().SetTitle("Mass (MeV/c^{2})")
    frame.Draw()
    frame.SetMinimum(0)
    canvas.SaveAs(f"fulldata_b7_down_time_binned_dpi_fit_bin_{i+1}.png")
    frame.Draw()
    canvas.SetLogy()
    frame.SetMinimum(3000)
    canvas.SaveAs(f"fulldata_b7_down_time_binned_dpi_fit_bin_{i+1}_log.png")

    print(f"Yield for {datetime.fromtimestamp(time_min / 1e6)} - {datetime.fromtimestamp(time_max / 1e6)}: {signal_yield:.2f}\n\n")# ± {signal_error:.2f}")

print(f"All yields: {signal_yields_dpi_down_b7}")
print(f"Background yields: {background_yields_dpi_down_b7}")


# # Time binned mass, $B^+ \rightarrow J/\psi K^+$

# ## Up

# In[ ]:


# import the data
tree_mjpsik_up = TChain("ST-b2cc")
tree_mjpsik_up.Add("more_more_data/combined/0028*_1.highstats-Small-B2CC-UP.root") # path to data

if tree_mjpsik_up.GetEntries() == 0:
    print("TChain is empty. Please check the file paths.")
    exit()


# In[ ]:


# create a RooRealVar for the observable to fit
x = r.RooRealVar("Bp_M", "Bp_M", 5160, 5400) #5200

# this to have time in the dataset for later
min_time_jpsik_up = tree_mjpsik_up.GetMinimum("GPSTIME")
max_time_jpsik_up = tree_mjpsik_up.GetMaximum("GPSTIME")
print(min_time_jpsik_up, max_time_jpsik_up)
print(datetime.fromtimestamp(min_time_jpsik_up / 1e6), datetime.fromtimestamp(max_time_jpsik_up / 1e6))
time = r.RooRealVar("GPSTIME", "GPSTIME", min_time_jpsik_up, max_time_jpsik_up)
block = r.RooRealVar("block", "block", 1, 8)

# convert the data to RooDataSet
data_jpsik_up = r.RooDataSet("data", "dataset from tree", r.RooArgSet(x, time, block), RooFit.Import(tree_mjpsik_up))


# In[ ]:


data_jpsik_up_b5 = data_jpsik_up.reduce("block >= 5 && block < 6")
data_jpsik_up_b8 = data_jpsik_up.reduce("block >= 8 && block < 9")


# In[ ]:


print(data_jpsik_up_b5.numEntries())
print(data_jpsik_up_b8.numEntries())


# ### Plots

# In[ ]:


c = r.TCanvas()

tree_mjpsik_up.Draw("block >> histblock(8, 1, 9)")

c.Draw()


# In[ ]:


c = r.TCanvas()

tree_mjpsik_up.Draw("Bp_M >> histmass(100, 5100, 5500)", "block >= 5 && block < 6", "logy")

c.Draw()


# In[ ]:


c = r.TCanvas()

tree_mjpsik_up.Draw("Bp_M >> histmass(100, 5160, 5400)")

hist = r.gPad.GetPrimitive("histmass")

hist.GetXaxis().SetTitle("Mass (MeV/c^{2})")
hist.GetYaxis().SetTitle("Events")
hist.SetTitle(f"B^{{+}} mass, B^{{+}} #rightarrow J/#psiK^{{+}}")

hist.SetStats(False)

r.gPad.Update()

c.Draw()


# ### Block 5

# In[ ]:


signal_yields_jpsik_up_b5 = np.empty(10, dtype=object)

background_yields_jpsik_up_b5 = np.empty(10, dtype=object)


# In[ ]:


for i in range(n_bins):
    time_min = time_bins_up_b5[i]
    time_max = time_bins_up_b5[i+1]

    # define the cut string for this time bin
    cut_str = f"GPSTIME >= {time_min} && GPSTIME < {time_max}"

    # create the binned dataset
    binned_data_jpsik_up = data_jpsik_up_b5.reduce(cut_str)

    print(f"Bin {i+1}: {time_min} - {time_max}, entries: {binned_data_jpsik_up.numEntries()}")

    # create the pdfs (2 crystal balls)
    mean = r.RooRealVar("mean", "mean", 5278.67, 5277, 5280)

    tau = r.RooRealVar("tau", "Decay constant", -0.000980333, -0.001, -0.0005)
    background = r.RooExponential("background", "Exponential background", x, tau)

    alpha_1 = r.RooRealVar("alpha", "alpha", 0.942242, 0.5, 2.5) # 1.10583
    n_1 = r.RooRealVar("n", "n", 7, 5, 10)
    cb_sigma_1 = r.RooRealVar("cb_sigma", "cb_sigma", 19.93, 15, 20)
    crystal_ball = r.RooCBShape("crystal_ball", "Crystal ball PDF", x, mean, cb_sigma_1, alpha_1, n_1)

    alpha_2 = r.RooRealVar("alpha2", "alpha2", -1.10518, -1.5, -0.7)
    n_2 = r.RooRealVar("n2", "n2", 50, 40, 60)
    cb_sigma_2 = r.RooRealVar("cb_sigma2", "cb_sigma2", 15, 10, 20)
    crystal_ball_2 = r.RooCBShape("crystal_ball2", "Crystal ball PDF 2", x, mean, cb_sigma_2, alpha_2, n_2)

    # fractions
    #frac_cb_2 = r.RooRealVar("frac_cb_2", "Fraction of crystal ball 2", 0.0537002, 0.0, 1.0)
    #frac_cb_1 = r.RooRealVar("frac_background", "Fraction of background", 0.199095, 0.0, 1.0)

    mean.setConstant(True)
    tau.setConstant(True)
    alpha_1.setConstant(False)
    n_1.setConstant(True)
    cb_sigma_1.setConstant(True)
    alpha_2.setConstant(True)
    n_2.setConstant(True)
    cb_sigma_2.setConstant(False)
    #frac_cb_2.setConstant(True)
    #frac_cb_1.setConstant(True)

    # use yields instead of fractions
    ntot = r.RooRealVar("ntot", "total yield",  binned_data_jpsik_up.numEntries(), 0,  binned_data_jpsik_up.numEntries() * 10)

#    ncb1 = r.RooRealVar("ncb1", "yield cb1", binned_data_dpi_up.numEntries() / 10, binned_data_dpi_up.numEntries() / 40, binned_data_dpi_up.numEntries())
#    ncb2 = r.RooRealVar("ncb2", "yield cb2", binned_data_dpi_up.numEntries() / 10, binned_data_dpi_up.numEntries() / 40, binned_data_dpi_up.numEntries())
    nsig = r.RooRealVar("nsig", "yield sig", binned_data_jpsik_up.numEntries() * 0.075, binned_data_jpsik_up.numEntries() * 0.06, binned_data_jpsik_up.numEntries() * 0.09)
    nbkg = r.RooRealVar("nbkg", "yield bkg", binned_data_jpsik_up.numEntries() * 0.915, binned_data_jpsik_up.numEntries() * 0.9, binned_data_jpsik_up.numEntries() * 0.93)
    f_cb1 = r.RooRealVar("f_cb1", "fraction of cb1", 0.6)

    # build the model and fit the data

    cb_sum = r.RooAddPdf("cb_sum", "signal pdf", r.RooArgList(crystal_ball, crystal_ball_2), r.RooArgList(f_cb1))

    combined_pdf_jpsik_up_timebinned = r.RooAddPdf("combined_pdf", "2x crystal ball + background",
                               r.RooArgList(cb_sum, background),
                               r.RooArgList(nsig, nbkg))

    # fit the model
    fit_result = combined_pdf_jpsik_up_timebinned.fitTo(binned_data_jpsik_up, RooFit.Save(),
    RooFit.Extended(True),         # use extended ML
    Strategy=2#,                   # more thorough minimization
    #RooFit.Hesse(True),           # always compute the Hesse matrix
    #RooFit.Minos(True)            # asymmetric uncertainties
    )

    print(x, mean, cb_sigma_1, alpha_1, n_1, cb_sigma_2, alpha_2, n_2, tau, ntot, nsig, nbkg)
    print(f"Background ratio: {nbkg.getVal() / ntot.getVal()}")
#    print(f"CB1/nsig ratio: {ncb1.getVal() / (ncb1.getVal() + ncb2.getVal())}")
    print(f"Status: {fit_result.status()}, CovQual: {fit_result.covQual()}")
    fit_result.Print()

    # extract the yield
    #nsig = r.RooAddition("nsig", "ncb1 + ncb2", r.RooArgList(ncb1, ncb2))

    #ncb1_err = ufloat(ncb1.getVal(), ncb1.getError())
    #ncb2_err = ufloat(ncb2.getVal(), ncb2.getError())

    #signal_yield = ncb1_err + ncb2_err

    signal_yield = ufloat(nsig.getVal(), nsig.getError())
    background_yield = ufloat(nbkg.getVal(), nbkg.getError())

    signal_yields_jpsik_up_b5[i] = signal_yield
    background_yields_jpsik_up_b5[i] = background_yield

    # plotting
    #frame = x.frame(RooFit.Title(f"Block 5 up JpsiK, time bin {i+1}: {datetime.fromtimestamp(time_min / 1e6)} - {datetime.fromtimestamp(time_max / 1e6)}"))
    frame = x.frame(RooFit.Title(f"B^{{+}} mass fit, B^{{+}} #rightarrow J/#psiK^{{+}}, block 5 (up), time bin {i+1}: {datetime.fromtimestamp(time_min / 1e6).date()} - {datetime.fromtimestamp(time_max / 1e6).date()}"))

    binned_data_jpsik_up.plotOn(frame)
    combined_pdf_jpsik_up_timebinned.plotOn(frame)
    combined_pdf_jpsik_up_timebinned.plotOn(frame, RooFit.Components("crystal_ball"), RooFit.LineStyle(r.kDashed), RooFit.LineColor(r.kBlue))
    combined_pdf_jpsik_up_timebinned.plotOn(frame, RooFit.Components("crystal_ball2"), RooFit.LineStyle(r.kDashed), RooFit.LineColor(r.kRed))
    combined_pdf_jpsik_up_timebinned.plotOn(frame, RooFit.Components("background"), RooFit.LineStyle(3), RooFit.LineColor(r.kGreen))

    #frame.GetYaxis().SetNdivisions(40)

    canvas = r.TCanvas(f"c_bin_{i+1}", f"Fit for time bin {i+1}", 800, 600)
    frame.GetXaxis().SetTitle("Mass (MeV/c^{2})")
    frame.Draw()
    frame.SetMinimum(0)
    canvas.SaveAs(f"fulldata_b5_up_time_binned_jpsik_fit_bin_{i+1}.png")
    frame.Draw()
    canvas.SetLogy()
    #canvas.SetTicky(1)
    frame.GetYaxis().SetMoreLogLabels(True)
    frame.GetYaxis().SetNdivisions(100520, True)

    latex = r.TLatex()
    latex.SetTextAngle(90)
    latex.SetTextSize(0.035)
    latex.SetTextFont(42)
    latex.SetTextAlign(22)  # center alignment
    latex.DrawLatexNDC(0.04, 0.55, "Events / ( 2.4 )")  # (x, y) in NDC
    frame.GetYaxis().SetTitle("")  # hide original title

    frame.GetXaxis().SetTitle("Mass (MeV/c^{2})")

    frame.SetMinimum(9999)
    canvas.SaveAs(f"fulldata_b5_up_time_binned_jpsik_fit_bin_{i+1}_log.png", "png 1000")

    print(f"Yield for {datetime.fromtimestamp(time_min / 1e6)} - {datetime.fromtimestamp(time_max / 1e6)}: {signal_yield:.2f}\n\n")# ± {signal_error:.2f}")

print(f"All yields: {signal_yields_jpsik_up_b5}")
print(f"Background yields: {background_yields_jpsik_up_b5}")


# ### Block 8

# In[ ]:


signal_yields_jpsik_up_b8 = np.empty(10, dtype=object)

background_yields_jpsik_up_b8 = np.empty(10, dtype=object)


# In[ ]:


for i in range(n_bins):
    time_min = time_bins_up_b8[i]
    time_max = time_bins_up_b8[i+1]

    # define the cut string for this time bin
    cut_str = f"GPSTIME >= {time_min} && GPSTIME < {time_max}"

    # create the binned dataset
    binned_data_jpsik_up = data_jpsik_up_b8.reduce(cut_str)

    print(f"Bin {i+1}: {time_min} - {time_max}, entries: {binned_data_jpsik_up.numEntries()}")

    # create the pdfs (2 crystal balls)
    mean = r.RooRealVar("mean", "mean", 5278.67, 5277, 5280)

    tau = r.RooRealVar("tau", "Decay constant", -0.000980333, -0.001, -0.0005)
    background = r.RooExponential("background", "Exponential background", x, tau)

    alpha_1 = r.RooRealVar("alpha", "alpha", 0.942242, 0.5, 2.5) # 1.10583
    n_1 = r.RooRealVar("n", "n", 7, 5, 10)
    cb_sigma_1 = r.RooRealVar("cb_sigma", "cb_sigma", 19.93, 15, 20)
    crystal_ball = r.RooCBShape("crystal_ball", "Crystal ball PDF", x, mean, cb_sigma_1, alpha_1, n_1)

    alpha_2 = r.RooRealVar("alpha2", "alpha2", -1.10518, -1.5, -0.7)
    n_2 = r.RooRealVar("n2", "n2", 50, 40, 60)
    cb_sigma_2 = r.RooRealVar("cb_sigma2", "cb_sigma2", 15, 10, 20)
    crystal_ball_2 = r.RooCBShape("crystal_ball2", "Crystal ball PDF 2", x, mean, cb_sigma_2, alpha_2, n_2)

    # fractions
    #frac_cb_2 = r.RooRealVar("frac_cb_2", "Fraction of crystal ball 2", 0.0537002, 0.0, 1.0)
    #frac_cb_1 = r.RooRealVar("frac_background", "Fraction of background", 0.199095, 0.0, 1.0)

    mean.setConstant(True)
    tau.setConstant(True)
    alpha_1.setConstant(False)
    n_1.setConstant(True)
    cb_sigma_1.setConstant(True)
    alpha_2.setConstant(True)
    n_2.setConstant(True)
    cb_sigma_2.setConstant(False)
    #frac_cb_2.setConstant(True)
    #frac_cb_1.setConstant(True)

    # use yields instead of fractions
    ntot = r.RooRealVar("ntot", "total yield",  binned_data_jpsik_up.numEntries(), 0,  binned_data_jpsik_up.numEntries() * 10)

#    ncb1 = r.RooRealVar("ncb1", "yield cb1", binned_data_dpi_up.numEntries() / 10, binned_data_dpi_up.numEntries() / 40, binned_data_dpi_up.numEntries())
#    ncb2 = r.RooRealVar("ncb2", "yield cb2", binned_data_dpi_up.numEntries() / 10, binned_data_dpi_up.numEntries() / 40, binned_data_dpi_up.numEntries())
    nsig = r.RooRealVar("nsig", "yield sig", binned_data_jpsik_up.numEntries() * 0.045, binned_data_jpsik_up.numEntries() * 0.03, binned_data_jpsik_up.numEntries() * 0.08)
    nbkg = r.RooRealVar("nbkg", "yield bkg", binned_data_jpsik_up.numEntries() * 0.945, binned_data_jpsik_up.numEntries() * 0.93, binned_data_jpsik_up.numEntries() * 0.97)
    f_cb1 = r.RooRealVar("f_cb1", "fraction of cb1", 0.6)

    # build the model and fit the data

    cb_sum = r.RooAddPdf("cb_sum", "signal pdf", r.RooArgList(crystal_ball, crystal_ball_2), r.RooArgList(f_cb1))

    combined_pdf_jpsik_up_timebinned = r.RooAddPdf("combined_pdf", "2x crystal ball + background",
                               r.RooArgList(cb_sum, background),
                               r.RooArgList(nsig, nbkg))

    # fit the model
    fit_result = combined_pdf_jpsik_up_timebinned.fitTo(binned_data_jpsik_up, RooFit.Save(),
    RooFit.Extended(True),         # use extended ML
    Strategy=2#,                   # more thorough minimization
    #RooFit.Hesse(True),           # always compute the Hesse matrix
    #RooFit.Minos(True)            # asymmetric uncertainties
    )

    print(x, mean, cb_sigma_1, alpha_1, n_1, cb_sigma_2, alpha_2, n_2, tau, ntot, nsig, nbkg)
    print(f"Background ratio: {nbkg.getVal() / ntot.getVal()}")
#    print(f"CB1/nsig ratio: {ncb1.getVal() / (ncb1.getVal() + ncb2.getVal())}")
    print(f"Status: {fit_result.status()}, CovQual: {fit_result.covQual()}")
    fit_result.Print()

    # extract the yield
    #nsig = r.RooAddition("nsig", "ncb1 + ncb2", r.RooArgList(ncb1, ncb2))

    #ncb1_err = ufloat(ncb1.getVal(), ncb1.getError())
    #ncb2_err = ufloat(ncb2.getVal(), ncb2.getError())

    #signal_yield = ncb1_err + ncb2_err

    signal_yield = ufloat(nsig.getVal(), nsig.getError())

    #signal_yield = nsig#.getVal()
    #signal_error = nsig.getError() # error propagation to be done

    signal_yield = ufloat(nsig.getVal(), nsig.getError())
    background_yield = ufloat(nbkg.getVal(), nbkg.getError())

    signal_yields_jpsik_up_b8[i] = signal_yield
    background_yields_jpsik_up_b8[i] = background_yield

    # plotting
    frame = x.frame(RooFit.Title(f"B^{{+}} mass fit, B^{{+}} #rightarrow J/#psiK^{{+}}, block 8 (up), time bin {i+1}: {datetime.fromtimestamp(time_min / 1e6).date()} - {datetime.fromtimestamp(time_max / 1e6).date()}"))
    binned_data_jpsik_up.plotOn(frame)
    combined_pdf_jpsik_up_timebinned.plotOn(frame)
    combined_pdf_jpsik_up_timebinned.plotOn(frame, RooFit.Components("crystal_ball"), RooFit.LineStyle(r.kDashed), RooFit.LineColor(r.kBlue))
    combined_pdf_jpsik_up_timebinned.plotOn(frame, RooFit.Components("crystal_ball2"), RooFit.LineStyle(r.kDashed), RooFit.LineColor(r.kRed))
    combined_pdf_jpsik_up_timebinned.plotOn(frame, RooFit.Components("background"), RooFit.LineStyle(3), RooFit.LineColor(r.kGreen))

    canvas = r.TCanvas(f"c_bin_{i+1}", f"Fit for time bin {i+1}", 800, 600)
    frame.GetXaxis().SetTitle("Mass (MeV/c^{2})")
    frame.Draw()
    frame.SetMinimum(0)
    canvas.SaveAs(f"fulldata_b8_up_time_binned_jpsik_fit_bin_{i+1}.png")
    frame.Draw()
    canvas.SetLogy()
    frame.SetMinimum(7000)
    canvas.SaveAs(f"fulldata_b8_up_time_binned_jpsik_fit_bin_{i+1}_log.png")

    print(f"Yield for {datetime.fromtimestamp(time_min / 1e6)} - {datetime.fromtimestamp(time_max / 1e6)}: {signal_yield:.2f}\n\n")# ± {signal_error:.2f}")

print(f"All yields: {signal_yields_jpsik_up_b8}")
print(f"Background yields: {background_yields_jpsik_up_b8}")


# ### JpsiK new fit (block 5)

# In[ ]:


#data_jpsik_up_b5 = data_jpsik_up.reduce("block >= 5 && block < 6")


# In[ ]:


# create the pdfs (2 crystal balls)
mean = r.RooRealVar("mean", "mean", 5278.67, 5277, 5280)

tau = r.RooRealVar("tau", "Decay constant", -0.000980333, -0.001, -0.0005)
background = r.RooExponential("background", "Exponential background", x, tau)

alpha_1 = r.RooRealVar("alpha", "alpha", 0.942242, 0.7, 1.5) # 1.10583
n_1 = r.RooRealVar("n", "n", 7, 5, 10)
cb_sigma_1 = r.RooRealVar("cb_sigma", "cb_sigma", 19.93, 15, 20)
crystal_ball = r.RooCBShape("crystal_ball", "Crystal ball PDF", x, mean, cb_sigma_1, alpha_1, n_1)

alpha_2 = r.RooRealVar("alpha2", "alpha2", -1.10518, -1.5, -0.7)
n_2 = r.RooRealVar("n2", "n2", 50, 40, 60)
cb_sigma_2 = r.RooRealVar("cb_sigma2", "cb_sigma2", 15, 10, 20)
crystal_ball_2 = r.RooCBShape("crystal_ball2", "Crystal ball PDF 2", x, mean, cb_sigma_2, alpha_2, n_2)

# fractions
#frac_cb_2 = r.RooRealVar("frac_cb_2", "Fraction of crystal ball 2", 0.0537002, 0.0, 1.0)
#frac_cb_1 = r.RooRealVar("frac_background", "Fraction of background", 0.199095, 0.0, 1.0)

mean.setConstant(True)
tau.setConstant(True)
alpha_1.setConstant(True)
n_1.setConstant(False)
cb_sigma_1.setConstant(True)
alpha_2.setConstant(True)
n_2.setConstant(True)
cb_sigma_2.setConstant(True)
#frac_cb_2.setConstant(True)
#frac_cb_1.setConstant(True)

# use yields instead of fractions
ntot = r.RooRealVar("ntot", "total yield",  data_jpsik_up_b5.numEntries(), 0,  data_jpsik_up_b5.numEntries() * 10)

ncb1 = r.RooRealVar("ncb1", "yield cb1", data_jpsik_up_b5.numEntries() / 50, data_jpsik_up_b5.numEntries() / 100, data_jpsik_up_b5.numEntries() / 10)
ncb2 = r.RooRealVar("ncb2", "yield cb2", data_jpsik_up_b5.numEntries() / 50, data_jpsik_up_b5.numEntries() / 100, data_jpsik_up_b5.numEntries() / 10)
nbkg = r.RooRealVar("nbkg", "yield bkg", data_jpsik_up_b5.numEntries() * 0.94, data_jpsik_up_b5.numEntries() * 0.93, data_jpsik_up_b5.numEntries() * 0.95)

# build the model and fit the data
combined_pdf_jpsik_up_timebinned = r.RooAddPdf("combined_pdf", "2x crystal ball + background",
                           r.RooArgList(crystal_ball_2, crystal_ball, background),
                           r.RooArgList(ncb2, ncb1, nbkg))

# fit the model
fit_result = combined_pdf_jpsik_up_timebinned.fitTo(data_jpsik_up_b5, RooFit.Save()#,
#r.RooFit.Extended(True),       # use extended ML
#r.RooFit.Strategy(2),          # more thorough minimization
#r.RooFit.Hesse(True),          # always compute the Hesse matrix
#r.RooFit.Minos(True)           # skip slow MINOS
)

print(x, mean, cb_sigma_1, alpha_1, n_1, cb_sigma_2, alpha_2, n_2, tau, ntot, ncb2, ncb1, nbkg)
print(f"Background ratio: {nbkg.getVal() / ntot.getVal()}")
print(f"Status: {fit_result.status()}, CovQual: {fit_result.covQual()}")
fit_result.Print()

# extract the yield
#nsig = r.RooAddition("nsig", "ncb1 + ncb2", r.RooArgList(ncb1, ncb2))

#ncb1_err = ufloat(ncb1.getVal(), ncb1.getError())
#ncb2_err = ufloat(ncb2.getVal(), ncb2.getError())

#signal_yield = ncb1_err + ncb2_err

#signal_yield = nsig#.getVal()
#signal_error = nsig.getError() # error propagation to be done

#signal_yields_jpsik_up[i] = signal_yield

# plotting
frame = x.frame(RooFit.Title("2CB JpsiK mass, block 5"))
data_jpsik_up_b5.plotOn(frame)
combined_pdf_jpsik_up_timebinned.plotOn(frame)
combined_pdf_jpsik_up_timebinned.plotOn(frame, RooFit.Components("crystal_ball"), RooFit.LineStyle(r.kDashed), RooFit.LineColor(r.kBlue))
combined_pdf_jpsik_up_timebinned.plotOn(frame, RooFit.Components("crystal_ball2"), RooFit.LineStyle(r.kDashed), RooFit.LineColor(r.kRed))
combined_pdf_jpsik_up_timebinned.plotOn(frame, RooFit.Components("background"), RooFit.LineStyle(3), RooFit.LineColor(r.kGreen))

canvas = r.TCanvas("canvas", "2x crystal ball fit", 800, 600)
frame.Draw()
canvas.SetLogy()
frame.SetMinimum(150000)
canvas.Draw()


# ## Down

# In[ ]:


# import the data
tree_mjpsik_down = TChain("ST-b2cc")
tree_mjpsik_down.Add("more_more_data/combined/0028*_1.highstats-Small-B2CC-DOWN.root") # path to data

if tree_mjpsik_down.GetEntries() == 0:
    print("TChain is empty. Please check the file paths.")
    exit()


# In[ ]:


# create a RooRealVar for the observable to fit
x = r.RooRealVar("Bp_M", "Bp_M", 5160, 5400) #5200

# this to have time in the dataset for later
min_time_jpsik_down = tree_mjpsik_down.GetMinimum("GPSTIME")
max_time_jpsik_down = tree_mjpsik_down.GetMaximum("GPSTIME")
print(min_time_jpsik_down, max_time_jpsik_down)
print(datetime.fromtimestamp(min_time_jpsik_down / 1e6), datetime.fromtimestamp(max_time_jpsik_down / 1e6))
time = r.RooRealVar("GPSTIME", "GPSTIME", min_time_jpsik_down, max_time_jpsik_down)
block = r.RooRealVar("block", "block", 1, 8)

# convert the data to RooDataSet
data_jpsik_down = r.RooDataSet("data", "dataset from tree", r.RooArgSet(x, time, block), RooFit.Import(tree_mjpsik_down))


# In[ ]:


data_jpsik_down_b6 = data_jpsik_down.reduce("block >= 6 && block < 7")
data_jpsik_down_b7 = data_jpsik_down.reduce("block >= 7 && block < 8")


# In[ ]:


print(data_jpsik_down_b6.numEntries())
print(data_jpsik_down_b7.numEntries())


# ### Block 6

# In[ ]:


signal_yields_jpsik_down_b6 = np.empty(10, dtype=object)

background_yields_jpsik_down_b6 = np.empty(10, dtype=object)


# In[ ]:


for i in range(n_bins):
    time_min = time_bins_down_b6[i]
    time_max = time_bins_down_b6[i+1]

    # define the cut string for this time bin
    cut_str = f"GPSTIME >= {time_min} && GPSTIME < {time_max}"

    # create the binned dataset
    binned_data_jpsik_down = data_jpsik_down_b6.reduce(cut_str)

    print(f"Bin {i+1}: {time_min} - {time_max}, entries: {binned_data_jpsik_down.numEntries()}")

    # create the pdfs (2 crystal balls)
    mean = r.RooRealVar("mean", "mean", 5278.67, 5277, 5280)

    tau = r.RooRealVar("tau", "Decay constant", -0.000980333, -0.001, -0.0005)
    background = r.RooExponential("background", "Exponential background", x, tau)

    alpha_1 = r.RooRealVar("alpha", "alpha", 0.942242, 0.7, 1.5) # 1.10583
    n_1 = r.RooRealVar("n", "n", 7, 5, 10)
    cb_sigma_1 = r.RooRealVar("cb_sigma", "cb_sigma", 19.93, 17, 22)
    crystal_ball = r.RooCBShape("crystal_ball", "Crystal ball PDF", x, mean, cb_sigma_1, alpha_1, n_1)

    alpha_2 = r.RooRealVar("alpha2", "alpha2", -1.10518, -1.5, -0.7)
    n_2 = r.RooRealVar("n2", "n2", 50, 40, 60)
    cb_sigma_2 = r.RooRealVar("cb_sigma2", "cb_sigma2", 15.8, 10, 20)
    crystal_ball_2 = r.RooCBShape("crystal_ball2", "Crystal ball PDF 2", x, mean, cb_sigma_2, alpha_2, n_2)

    # fractions
    #frac_cb_2 = r.RooRealVar("frac_cb_2", "Fraction of crystal ball 2", 0.0537002, 0.0, 1.0)
    #frac_cb_1 = r.RooRealVar("frac_background", "Fraction of background", 0.199095, 0.0, 1.0)

    mean.setConstant(True)
    tau.setConstant(True)
    alpha_1.setConstant(False)
    n_1.setConstant(True)
    cb_sigma_1.setConstant(True)
    alpha_2.setConstant(True)
    n_2.setConstant(True)
    cb_sigma_2.setConstant(False)
    #frac_cb_2.setConstant(True)
    #frac_cb_1.setConstant(True)

    # use yields instead of fractions
    ntot = r.RooRealVar("ntot", "total yield",  binned_data_jpsik_down.numEntries(), 0,  binned_data_jpsik_down.numEntries() * 10)

    nsig = r.RooRealVar("nsig", "yield sig", binned_data_jpsik_down.numEntries() * 0.075, binned_data_jpsik_down.numEntries() * 0.06, binned_data_jpsik_down.numEntries() * 0.09)
    nbkg = r.RooRealVar("nbkg", "yield bkg", binned_data_jpsik_down.numEntries() * 0.915, binned_data_jpsik_down.numEntries() * 0.9, binned_data_jpsik_down.numEntries() * 0.93)
    f_cb1 = r.RooRealVar("f_cb1", "fraction of cb1", 0.6)

    # build the model and fit the data

    cb_sum = r.RooAddPdf("cb_sum", "signal pdf", r.RooArgList(crystal_ball, crystal_ball_2), r.RooArgList(f_cb1))

    combined_pdf_jpsik_down_timebinned = r.RooAddPdf("combined_pdf", "2x crystal ball + background",
                               r.RooArgList(cb_sum, background),
                               r.RooArgList(nsig, nbkg))

    # fit the model
    fit_result = combined_pdf_jpsik_down_timebinned.fitTo(binned_data_jpsik_down, RooFit.Save(),
    RooFit.Extended(True),         # use extended ML
    Strategy=2#,                   # more thorough minimization
    #RooFit.Hesse(True),           # always compute the Hesse matrix
    #RooFit.Minos(True)            # asymmetric uncertainties
    )

    print(x, mean, cb_sigma_1, alpha_1, n_1, cb_sigma_2, alpha_2, n_2, tau, ntot, nsig, nbkg)
    print(f"Background ratio: {nbkg.getVal() / ntot.getVal()}")
#    print(f"CB1/nsig ratio: {ncb1.getVal() / (ncb1.getVal() + ncb2.getVal())}")
    print(f"Status: {fit_result.status()}, CovQual: {fit_result.covQual()}")
    fit_result.Print()

    signal_yield = ufloat(nsig.getVal(), nsig.getError())
    background_yield = ufloat(nbkg.getVal(), nbkg.getError())

    signal_yields_jpsik_down_b6[i] = signal_yield
    background_yields_jpsik_down_b6[i] = background_yield

    # plotting
    frame = x.frame(RooFit.Title(f"B^{{+}} mass fit, B^{{+}} #rightarrow J/#psiK^{{+}}, block 6 (down)), time bin {i+1}: {datetime.fromtimestamp(time_min / 1e6).date()} - {datetime.fromtimestamp(time_max / 1e6).date()}"))
    binned_data_jpsik_down.plotOn(frame)
    combined_pdf_jpsik_down_timebinned.plotOn(frame)
    combined_pdf_jpsik_down_timebinned.plotOn(frame, RooFit.Components("crystal_ball"), RooFit.LineStyle(r.kDashed), RooFit.LineColor(r.kBlue))
    combined_pdf_jpsik_down_timebinned.plotOn(frame, RooFit.Components("crystal_ball2"), RooFit.LineStyle(r.kDashed), RooFit.LineColor(r.kRed))
    combined_pdf_jpsik_down_timebinned.plotOn(frame, RooFit.Components("background"), RooFit.LineStyle(3), RooFit.LineColor(r.kGreen))

    canvas = r.TCanvas(f"c_bin_{i+1}", f"Fit for time bin {i+1}", 800, 600)
    frame.GetXaxis().SetTitle("Mass (MeV/c^{2})")
    frame.Draw()
    frame.SetMinimum(0)
    canvas.SaveAs(f"fulldata_b6_down_time_binned_jpsik_fit_bin_{i+1}.png")
    frame.Draw()
    canvas.SetLogy()
    frame.GetYaxis().SetMoreLogLabels(True)
    frame.GetYaxis().SetNdivisions(100520, True)

    latex = r.TLatex()
    latex.SetTextAngle(90)
    latex.SetTextSize(0.035)
    latex.SetTextFont(42)
    latex.SetTextAlign(22)  # center alignment
    latex.DrawLatexNDC(0.04, 0.55, "Events / ( 2.4 )")  # (x, y) in NDC
    frame.GetYaxis().SetTitle("")  # hide original title

    frame.GetXaxis().SetTitle("Mass (MeV/c^{2})")

    frame.SetMinimum(9999)
    frame.SetMaximum(20001)
    canvas.SaveAs(f"fulldata_b6_down_time_binned_jpsik_fit_bin_{i+1}_log.png")

    print(f"Yield for {datetime.fromtimestamp(time_min / 1e6)} - {datetime.fromtimestamp(time_max / 1e6)}: {signal_yield:.2f}\n\n")# ± {signal_error:.2f}")

print(f"All yields: {signal_yields_jpsik_down_b6}")
print(f"Background yields: {background_yields_jpsik_down_b6}")


# ### Block 7

# In[ ]:


signal_yields_jpsik_down_b7 = np.empty(10, dtype=object)

background_yields_jpsik_down_b7 = np.empty(10, dtype=object)


# In[ ]:


for i in range(n_bins):
    time_min = time_bins_down_b7[i]
    time_max = time_bins_down_b7[i+1]

    # define the cut string for this time bin
    cut_str = f"GPSTIME >= {time_min} && GPSTIME < {time_max}"

    # create the binned dataset
    binned_data_jpsik_down = data_jpsik_down_b7.reduce(cut_str)

    print(f"Bin {i+1}: {time_min} - {time_max}, entries: {binned_data_jpsik_down.numEntries()}")

    # create the pdfs (2 crystal balls)
    mean = r.RooRealVar("mean", "mean", 5278.67, 5277, 5280)

    tau = r.RooRealVar("tau", "Decay constant", -0.000984097, -0.001, -0.0005)
    background = r.RooExponential("background", "Exponential background", x, tau)

    # fix misbehaving bin, could be wise to figure out why it was misbehaving instead
    if i == 8:
        alpha_1 = r.RooRealVar("alpha", "alpha", 1.47924, 0.5, 2.5)
        alpha_1.setConstant(True)

    else:
        alpha_1 = r.RooRealVar("alpha", "alpha", 1.20315, 0.5, 2.5) # 1.10583
        alpha_1.setConstant(False)

    n_1 = r.RooRealVar("n", "n", 3.62277, 3, 10)
    cb_sigma_1 = r.RooRealVar("cb_sigma", "cb_sigma", 19.6, 15, 25)
    crystal_ball = r.RooCBShape("crystal_ball", "Crystal ball PDF", x, mean, cb_sigma_1, alpha_1, n_1)

    alpha_2 = r.RooRealVar("alpha2", "alpha2", -1.13179, -1.5, -0.8)
    n_2 = r.RooRealVar("n2", "n2", 72.0858, 40, 100)
    cb_sigma_2 = r.RooRealVar("cb_sigma2", "cb_sigma2", 16.5, 14, 19)
    crystal_ball_2 = r.RooCBShape("crystal_ball2", "Crystal ball PDF 2", x, mean, cb_sigma_2, alpha_2, n_2)

    # fractions
    #frac_cb_2 = r.RooRealVar("frac_cb_2", "Fraction of crystal ball 2", 0.0537002, 0.0, 1.0)
    #frac_cb_1 = r.RooRealVar("frac_background", "Fraction of background", 0.199095, 0.0, 1.0)

    mean.setConstant(True)
    tau.setConstant(True)
    #alpha_1.setConstant(False)
    n_1.setConstant(True)
    cb_sigma_1.setConstant(False)
    alpha_2.setConstant(True)
    n_2.setConstant(True)
    cb_sigma_2.setConstant(True)
    #frac_cb_2.setConstant(True)
    #frac_cb_1.setConstant(True)

    # use yields instead of fractions
    ntot = r.RooRealVar("ntot", "total yield",  binned_data_jpsik_down.numEntries(), 0,  binned_data_jpsik_down.numEntries() * 10)

    nsig = r.RooRealVar("nsig", "yield sig", binned_data_jpsik_down.numEntries() * 0.075, binned_data_jpsik_down.numEntries() * 0.03, binned_data_jpsik_down.numEntries() * 0.09)
    nbkg = r.RooRealVar("nbkg", "yield bkg", binned_data_jpsik_down.numEntries() * 0.95, binned_data_jpsik_down.numEntries() * 0.915, binned_data_jpsik_down.numEntries() * 0.97)
    f_cb1 = r.RooRealVar("f_cb1", "fraction of cb1", 0.6)

    # build the model and fit the data

    cb_sum = r.RooAddPdf("cb_sum", "signal pdf", r.RooArgList(crystal_ball, crystal_ball_2), r.RooArgList(f_cb1))

    combined_pdf_jpsik_down_timebinned = r.RooAddPdf("combined_pdf", "2x crystal ball + background",
                               r.RooArgList(cb_sum, background),
                               r.RooArgList(nsig, nbkg))

   # alpha.setError(0.01)

    # fit the model
    fit_result = combined_pdf_jpsik_down_timebinned.fitTo(binned_data_jpsik_down, RooFit.Save(),
    RooFit.Extended(True),         # use extended ML
    RooFit.Strategy(2),            # more thorough minimization
    #RooFit.Hesse(True),           # always compute the Hesse matrix
    #RooFit.Minos(True)            # asymmetric uncertainties
    RooFit.Minimizer("Minuit2"),
    RooFit.MaxCalls(2000000000)
    )

    print(x, mean, cb_sigma_1, alpha_1, n_1, cb_sigma_2, alpha_2, n_2, tau, ntot, nsig, nbkg)
    print(f"Background ratio: {nbkg.getVal() / ntot.getVal()}")
#    print(f"CB1/nsig ratio: {ncb1.getVal() / (ncb1.getVal() + ncb2.getVal())}")
    print(f"Status: {fit_result.status()}, CovQual: {fit_result.covQual()}")
    fit_result.Print()

    signal_yield = ufloat(nsig.getVal(), nsig.getError())
    background_yield = ufloat(nbkg.getVal(), nbkg.getError())

    signal_yields_jpsik_down_b7[i] = signal_yield
    background_yields_jpsik_down_b7[i] = background_yield

    # plotting
    frame = x.frame(RooFit.Title(f"B^{{+}} mass fit, B^{{+}} #rightarrow J/#psiK^{{+}}, block 7 (down), time bin {i+1}: {datetime.fromtimestamp(time_min / 1e6).date()} - {datetime.fromtimestamp(time_max / 1e6).date()}"))
    binned_data_jpsik_down.plotOn(frame)
    combined_pdf_jpsik_down_timebinned.plotOn(frame)
    combined_pdf_jpsik_down_timebinned.plotOn(frame, RooFit.Components("crystal_ball"), RooFit.LineStyle(r.kDashed), RooFit.LineColor(r.kBlue))
    combined_pdf_jpsik_down_timebinned.plotOn(frame, RooFit.Components("crystal_ball2"), RooFit.LineStyle(r.kDashed), RooFit.LineColor(r.kRed))
    combined_pdf_jpsik_down_timebinned.plotOn(frame, RooFit.Components("background"), RooFit.LineStyle(3), RooFit.LineColor(r.kGreen))

    canvas = r.TCanvas(f"c_bin_{i+1}", f"Fit for time bin {i+1}", 800, 600)
    frame.GetXaxis().SetTitle("Mass (MeV/c^{2})")
    frame.Draw()
    frame.SetMinimum(0)
    canvas.SaveAs(f"fulldata_b7_down_time_binned_jpsik_fit_bin_{i+1}.png")
    frame.Draw()
    canvas.SetLogy()
    frame.GetYaxis().SetMoreLogLabels(True)
    frame.GetYaxis().SetNdivisions(100520, True)

    latex = r.TLatex()
    latex.SetTextAngle(90)
    latex.SetTextSize(0.035)
    latex.SetTextFont(42)
    latex.SetTextAlign(22)  # center alignment
    latex.DrawLatexNDC(0.04, 0.55, "Events / ( 2.4 )")  # (x, y) in NDC
    frame.GetYaxis().SetTitle("")  # hide original title

    frame.GetXaxis().SetTitle("Mass (MeV/c^{2})")

    frame.SetMinimum(9999)
    frame.SetMaximum(20001)
    canvas.SaveAs(f"fulldata_b7_down_time_binned_jpsik_fit_bin_{i+1}_log.png")

    print(f"Yield for {datetime.fromtimestamp(time_min / 1e6)} - {datetime.fromtimestamp(time_max / 1e6)}: {signal_yield:.2f}\n\n")# ± {signal_error:.2f}")

print(f"All yields: {signal_yields_jpsik_down_b7}")
print(f"Background yields: {background_yields_jpsik_down_b7}")


# # Ratios

# In[ ]:


signal_yields_dpi = [signal_yields_dpi_up_b5, signal_yields_dpi_down_b6, signal_yields_dpi_down_b7, signal_yields_dpi_up_b8]
signal_yields_jpsik = [signal_yields_jpsik_up_b5, signal_yields_jpsik_down_b6, signal_yields_jpsik_down_b7, signal_yields_jpsik_up_b8]
block_time_bins = [time_bins_up_b5, time_bins_down_b6, time_bins_down_b7, time_bins_up_b8]


# In[ ]:


print(signal_yields_dpi)
print(signal_yields_jpsik)


# In[ ]:


avg_ratios = np.empty(4, dtype=object)


# ## All blocks

# In[ ]:


for i in range(4):

    ratios = [a / b for a, b in zip(signal_yields_dpi[i], signal_yields_jpsik[i])]

    ratios_val = np.array([r.nominal_value for r in ratios])
    ratios_err = np.array([r.std_dev for r in ratios])

    print(signal_yields_dpi[i])
    print(signal_yields_jpsik[i])
    print(ratios)

    x_plot = np.arange(len(ratios_val))
    plt.errorbar(x_plot, ratios_val, yerr=ratios_err, fmt='o', capsize=5, label=r"$N_{\bar{D}^0\pi^+}$ / $N_{J/\psi K^+}$", color='black', ecolor='black')
    plt.xlabel("Time interval")
    plt.ylabel(r"$B^+ \rightarrow \bar{D}^0\pi^+$ / $B^+ \rightarrow J/\psi K^+$ ratio")
    if i == 0 or i == 3:
        plt.title(r"$N_{\bar{D}^0\pi^+}$ / $N_{J/\psi K^+}$" + f", block {i+5} (up)")
    else:
        plt.title(r"$N_{\bar{D}^0\pi^+}$ / $N_{J/\psi K^+}$" + f", block {i+5} (down)")
    ticks = [f"{datetime.fromtimestamp(block_time_bins[i][k] / 1e6).date()} - \n{datetime.fromtimestamp(block_time_bins[i][k+1] / 1e6).date()}" for k in range(n_bins)]
    plt.xticks(x_plot, ticks, rotation=60, ha='right', fontsize=8)
    #plt.ylim(0.95, 1.2)
    plt.legend()
    plt.grid(True)
    plt.show()

    # fit and chi^2

    constant_function = lambda x, a: a

    popt, pcov = curve_fit(constant_function, x_plot, ratios_val, sigma=ratios_err)#, absolute_sigma=True)

    y_fit = constant_function(x_plot, *popt)
    y_unc = np.sqrt(pcov[0, 0])

    y_fit_err = ufloat(y_fit, y_unc)
    avg_ratios[i] = y_fit_err

    # calculate chi^2
    chi_squared = np.sum(((ratios_val - y_fit) / ratios_err) ** 2)

    chi_squared_reduced = chi_squared / 9

    print(f"Fitted constant: {y_fit_err}")
    print(f"Chi^2: {chi_squared}")
    print(f"Chi^2/dof: {chi_squared_reduced}")

    plt.errorbar(x_plot, ratios_val, yerr=ratios_err, fmt='o', capsize=5, label=r"$N_{\bar{D}^0\pi^+}$ / $N_{J/\psi K^+}$", color='black', ecolor='black')
    plt.axhline(y_fit, linestyle="dotted", color='red', label=rf"Fit at $y={y_fit_err.nominal_value:.3f} \pm {math.ceil(y_fit_err.std_dev * 10 ** 3) / (10 ** 3)}$")
    plt.xlabel("Time interval")
    plt.ylabel(r"$B^+ \rightarrow \bar{D}^0\pi^+$ / $B^+ \rightarrow J/\psi K^+$ ratio")
    if i == 0 or i == 3:
        plt.title(rf"$N_{{\bar{{D}}^0\pi^+}}$ / $N_{{J/\psi K^+}}$, block {i+5} (up); $\chi^2 = {chi_squared:.2f}$, $\chi^2/\text{{dof}} = {chi_squared_reduced:.2f}$")
    else:
        plt.title(rf"$N_{{\bar{{D}}^0\pi^+}}$ / $N_{{J/\psi K^+}}$, block {i+5} (down); $\chi^2 = {chi_squared:.2f}$, $\chi^2/\text{{dof}} = {chi_squared_reduced:.2f}$")

    ticks = [f"{datetime.fromtimestamp(block_time_bins[i][k] / 1e6).date()} - \n{datetime.fromtimestamp(block_time_bins[i][k+1] / 1e6).date()}" for k in range(n_bins)]
    plt.xticks(x_plot, ticks, rotation=60, ha='right', fontsize=8)
    #plt.ylim(0.95, 1.2)
    if i == 0 or i == 2:
        plt.legend(loc="upper right")
    else:
        plt.legend(loc="upper left")
    plt.grid(True)
    plt.show()


# ## Yield ratios per block

# In[ ]:


avg_ratios_val = np.array([r.nominal_value for r in avg_ratios])
avg_ratios_err = np.array([r.std_dev for r in avg_ratios])

print(avg_ratios)

x_plot = np.arange(len(avg_ratios_val))
plt.errorbar(x_plot, avg_ratios_val, yerr=avg_ratios_err, fmt='o', capsize=5, label=r"Average $N_{\bar{D}^0\pi^+}$ / $N_{J/\psi K^+}$", color='black', ecolor='black')
#plt.xlabel("Block")
plt.ylabel(r"Fitted $B^+ \rightarrow \bar{D}^0\pi^+$ / $B^+ \rightarrow J/\psi K^+$ ratio")
plt.title(r"Average $N_{\bar{D}^0\pi^+}$ / $N_{J/\psi K^+}$" + f" per block")
ticks = [f"Block {i+5};\n{datetime.fromtimestamp(block_time_bins[i][0] / 1e6).date()} -\n{datetime.fromtimestamp(block_time_bins[i][-1] / 1e6).date()}" for i in range(4)]
plt.xticks(x_plot, ticks, fontsize=8)#, rotation=60, ha='right', fontsize=8)
#plt.ylim(0.95, 1.2)
plt.legend()
plt.grid(True)
plt.show()

# fit and chi^2

constant_function = lambda x, a: a

popt, pcov = curve_fit(constant_function, x_plot, avg_ratios_val, sigma=avg_ratios_err)#, absolute_sigma=True)

y_fit = constant_function(x_plot, *popt)
y_unc = np.sqrt(pcov[0, 0])

y_fit_err = ufloat(y_fit, y_unc)

# calculate chi^2
chi_squared = np.sum(((avg_ratios_val - y_fit) / avg_ratios_err) ** 2)

chi_squared_reduced = chi_squared / 3

print(f"Fitted constant: {y_fit_err}")
print(f"Chi^2: {chi_squared}")
print(f"Chi^2/dof: {chi_squared_reduced}")

plt.errorbar(x_plot, avg_ratios_val, yerr=avg_ratios_err, fmt='o', capsize=5, label=r"Average $N_{\bar{D}^0\pi^+}$ / $N_{J/\psi K^+}$", color='black', ecolor='black')
plt.axhline(y_fit, linestyle="dotted", color='red', label=rf"Fit at $y={y_fit_err.nominal_value:.3f} \pm {math.ceil(y_fit_err.std_dev * 10 ** 3) / (10 ** 3)}$")
#plt.xlabel("Block")
plt.ylabel(r"Fitted $B^+ \rightarrow \bar{D}^0\pi^+$ / $B^+ \rightarrow J/\psi K^+$ ratio")
plt.title(rf"Average $N_{{\bar{{D}}^0\pi^+}}$ / $N_{{J/\psi K^+}}$ per block")#; $\chi^2 = {chi_squared:.2f}$, $\chi^2/\text{{dof}} = {chi_squared_reduced:.2f}$")
ticks = [f"Block {i+5};\n{datetime.fromtimestamp(block_time_bins[i][0] / 1e6).date()} -\n{datetime.fromtimestamp(block_time_bins[i][-1] / 1e6).date()}" for i in range(4)]
plt.xticks(x_plot, ticks, fontsize=8)#, rotation=30, ha='right', fontsize=8)
#plt.ylim(0.95, 1.2)
plt.legend()
plt.grid(True)
plt.show()


# # Variable plotting (incomplete)

# I removed the work on plotting variables as they were not necessary in the end when the data proved (preliminarily) stable; only the complete SNR plots are left below.

# ## SNR

# ### Function

# In[ ]:


def plot_snr(signal_yields, background_yields, time_bins, decay, block, polarity, print_values=False):

    snr = [a / b for a, b in zip(signal_yields, background_yields)]

    ratios_val = np.array([r.nominal_value for r in snr])
    ratios_err = np.array([r.std_dev for r in snr])

    if print_values:
        print(signal_yields_dpi_up_b5)
        print(background_yields_dpi_up_b5)
        print(snr_dpi_up_b5)

    x_plot = np.arange(len(ratios_val))
    plt.errorbar(x_plot, ratios_val, yerr=ratios_err, fmt='o', capsize=5, label="SNR", color='black', ecolor='black')
    plt.xlabel("Time interval")
    plt.ylabel("SNR")
    if decay == "dpi":
        plt.title(rf"Signal to background ratio, $B^+ \rightarrow \bar{{D}}^0\pi^+$, block {block} ({polarity})")
    if decay == "jpsik":
        plt.title(rf"Signal to background ratio, $B^+ \rightarrow J/\psi K^+$, block {block} ({polarity})")
    ticks = [f"{datetime.fromtimestamp(time_bins[i] / 1e6).date()} - \n{datetime.fromtimestamp(time_bins[i+1] / 1e6).date()}" for i in range(n_bins)]
    plt.xticks(x_plot, ticks, rotation=60, ha='right', fontsize=8)
    #plt.ylim(0.95, 1.2)
    plt.legend()
    plt.grid(True)
    plt.show()


# ### $B^+ \rightarrow \bar{D}^0\pi^+$

# #### Block 5 (up)

# In[ ]:


plot_snr(signal_yields_dpi_up_b5,
         background_yields_dpi_up_b5,
         time_bins_up_b5,
         "dpi",
         5,
         "up",
         print_values=False)


# #### Block 6 (down)

# In[ ]:


plot_snr(signal_yields_dpi_down_b6, background_yields_dpi_down_b6, time_bins_down_b6, "dpi", 6, "down", print_values=False)


# #### Block 7 (down)

# In[ ]:


plot_snr(signal_yields_dpi_down_b7, background_yields_dpi_down_b7, time_bins_down_b7, "dpi", 7, "down", print_values=False)


# #### Block 8 (up)

# In[ ]:


plot_snr(signal_yields_dpi_up_b8, background_yields_dpi_up_b8, time_bins_up_b8, "dpi", 8, "up", print_values=False)


# ### $B^+ \rightarrow J/\psi K^+$

# #### Block 5 (up)

# In[ ]:


plot_snr(signal_yields_jpsik_up_b5, background_yields_jpsik_up_b5, time_bins_up_b5, "jpsik", 5, "up", print_values=False)


# #### Block 6 (down)

# In[ ]:


plot_snr(signal_yields_jpsik_down_b6, background_yields_jpsik_down_b6, time_bins_down_b6, "jpsik", 6, "down", print_values=False)


# #### Block 7 (down)

# In[ ]:


plot_snr(signal_yields_jpsik_down_b7, background_yields_jpsik_down_b7, time_bins_down_b7, "jpsik", 7, "down", print_values=False)


# #### Block 8 (up)

# In[ ]:


plot_snr(signal_yields_jpsik_up_b8, background_yields_jpsik_up_b8, time_bins_up_b8, "jpsik", 8, "up", print_values=False)

