"""Generate the two mathematical Model analysis figures."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from analysis_math import ROOT, GAMMA_B, LOWER_SLOPE, derivative, root, checks

FIGURES=ROOT/"figures"
FIGURES.mkdir(exist_ok=True)

def format_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.9)
    ax.spines["bottom"].set_linewidth(0.9)
    ax.tick_params(axis="both",direction="out",width=0.9,length=4,labelsize=11)

def save(fig,name):
    fig.savefig(FIGURES/f"{name}.eps",format="eps",bbox_inches="tight",pad_inches=0.08)
    fig.savefig(FIGURES/f"{name}.png",dpi=220,bbox_inches="tight",pad_inches=0.08)

def main():
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":12,
        "mathtext.fontset":"dejavusans","axes.labelsize":13,
        "legend.fontsize":10.5,"ps.fonttype":42,"lines.linewidth":1.8})
    report=checks()
    assert report["NBC_capacity_recovery_difference"]<1e-13
    assert report["WT"]["V_b_error_V"]<1e-12
    assert report["KO"]["V_b_error_V"]<1e-12
    assert report["maximum_relative_carrier_symmetry_error"]<1e-13
    (ROOT/"numerical_checks.json").write_text(json.dumps(report,indent=2)+"\n")

    potentials=np.linspace(-8.0,0.0,1201)
    fig,ax=plt.subplots(figsize=(7.0,4.6))
    table=[potentials]
    for case,label,linestyle in [("WT","Control","-"),("KO","Ae4 knockout","--")]:
        vals=derivative(case,potentials)
        ps=root(case)
        xx=np.sort(np.r_[potentials,ps])
        idx=int(np.flatnonzero(xx==ps)[0])
        ax.plot(xx,derivative(case,xx),linestyle=linestyle,label=label,
                marker="o",markevery=[idx],markersize=5.2)
        table.append(vals)
    ax.plot(potentials,np.full_like(potentials,LOWER_SLOPE),linestyle=":",
            linewidth=1.5,label="Positive channel contribution")
    ax.set_xlabel(r"Basolateral membrane potential, $\psi_b$")
    ax.set_ylabel(r"$d\mathcal{H}/d\psi_b$")
    ax.set_xlim(-8.0,0.0)
    ax.set_ylim(LOWER_SLOPE-0.09,LOWER_SLOPE+GAMMA_B/2+0.18)
    ax.xaxis.set_major_locator(MultipleLocator(1.0))
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.legend(frameon=False,loc="upper left",handlelength=2.5)
    ax.text(0.985,0.95,r"$\tau=1,\quad u=1$",transform=ax.transAxes,ha="right",va="top",fontsize=11)
    ax.text(0.985,0.82,"Ionic compositions held fixed",transform=ax.transAxes,ha="right",va="top",fontsize=10)
    ax.text(0.60,0.12,r"$g_b+g_p g_a/(g_a+g_p)=12.34>0$",transform=ax.transAxes,
            ha="center",va="bottom",fontsize=12,bbox=dict(facecolor="white",edgecolor="none",pad=1.8))
    format_axes(ax)
    fig.subplots_adjust(left=.13,right=.985,bottom=.16,top=.97)
    save(fig,"Ae4_voltage_uniqueness")
    np.savetxt(ROOT/"voltage_uniqueness_data.csv",np.column_stack(table),delimiter=",",
               header="psi_b,Hprime_control,Hprime_knockout",comments="")
    plt.close(fig)

    fig,ax=plt.subplots(figsize=(7.0,4.6))
    curves=[]
    for ratio in [1.0,0.5,2.0,4.0]:
        lower=max(.25,ratio/3.0); upper=min(3.0,ratio/.25)
        a=np.linspace(lower,upper,700); b=ratio/a
        ax.plot(a,b,linestyle="-" if ratio==1 else "--",linewidth=2.4 if ratio==1 else 1.5)
        curves.extend(np.column_stack((np.full(len(a),ratio),a,b)).tolist())
    pts=np.array([[.5,2.],[1.,1.],[2.,.5]])
    ax.plot(pts[:,0],pts[:,1],linestyle="None",marker="o",markersize=6)
    ax.annotate(r"$(0.5,\,2)$",(.5,2.),xytext=(9,8),textcoords="offset points",fontsize=10.5)
    ax.annotate(r"Reference $(1,\,1)$",(1.,1.),xytext=(11,9),textcoords="offset points",fontsize=10.5)
    ax.annotate(r"$(2,\,0.5)$",(2.,.5),xytext=(8,12),textcoords="offset points",fontsize=10.5)
    ax.text(.985,.97,r"$r_{\mathrm{Na}}$ and $r_{\mathrm{K}}$ fixed",transform=ax.transAxes,ha="right",va="top",fontsize=11)
    ax.text(.985,.87,"Identical trajectories\nalong each curve",transform=ax.transAxes,ha="right",va="top",fontsize=10.5)
    ax.text(2.92,1/2.92+.04,r"$\Lambda_4/\Lambda_{4,\mathrm{ref}}=1$",ha="right",va="bottom",fontsize=11)
    ax.text(1.35,.5/1.35+.06,r"$0.5$",ha="left",va="bottom",fontsize=11)
    ax.text(2.84,2/2.84+.055,r"$2$",ha="right",va="bottom",fontsize=11)
    ax.text(2.84,4/2.84+.055,r"$4$",ha="right",va="bottom",fontsize=11)
    ax.set_xlim(.25,3.05); ax.set_ylim(.22,3.05)
    ax.set_xlabel(r"Carrier amount, $N_4/N_{4,\mathrm{ref}}$")
    ax.set_ylabel(r"Common rate scale, $k_0/k_{0,\mathrm{ref}}$")
    ax.xaxis.set_major_locator(MultipleLocator(.5)); ax.yaxis.set_major_locator(MultipleLocator(.5))
    format_axes(ax)
    fig.subplots_adjust(left=.13,right=.985,bottom=.16,top=.97)
    save(fig,"Ae4_parameter_equivalence")
    np.savetxt(ROOT/"parameter_equivalence_data.csv",np.asarray(curves),delimiter=",",
               header="Lambda4_ratio,carrier_amount_ratio,common_rate_ratio",comments="")
    plt.close(fig)

if __name__=="__main__":
    main()
