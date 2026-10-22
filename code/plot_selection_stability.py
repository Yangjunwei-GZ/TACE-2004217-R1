from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
def plot_stability(summary,destination):
    with matplotlib.rc_context(matplotlib.rcParamsDefault):
        plt.rcParams.update({'font.family':'Arial','font.size':9,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
        fig,axs=plt.subplots(1,2,figsize=(7.2,3.35),sharex=True,sharey=True)
        labels=['SZN','LAHGLE','Both selected']
        for ax,mode,panel in zip(axs,['point','bootstrap'],['a','b']):
         rows=[x for x in summary if x['mode']==mode]
         vals=[x['percent'] for x in rows]
         ax.barh([2,1,0],vals,height=.53,color=['#537B94','#7899AD','#B2B9BE'])
         ax.set_yticks([2,1,0],labels);ax.set_xlim(0,100);ax.set_xticks([0,25,50,75,100]);ax.set_xlabel('Selection frequency (%)')
         ax.set_title(f'{panel}   '+('Point: 50 training fits' if mode=='point' else 'Bootstrap: 2,500 training fits'),loc='left',fontsize=10,pad=14)
         for y,r in zip([2,1,0],rows):ax.text(r['percent']+2,y,f"{r['selected_fits']}/{r['total_fits']:,}\n({r['percent']:.2f}%)",va='center',fontsize=8)
         ax.grid(axis='x',color='#E5E8EB',linewidth=.6);ax.set_axisbelow(True)
        fig.subplots_adjust(left=.14,right=.98,bottom=.26,top=.84,wspace=.35)
        fig.text(.14,.085,'Selected: |coefficient| > 1e−8. Fits overlap; no independent-binomial intervals.',fontsize=8)
        for ext in ['svg','pdf','tiff','png']:fig.savefig(str(destination)+'.'+ext,dpi=600 if ext=='tiff' else 160,facecolor='white')
        plt.close(fig)
