# ===================================================================
# Created: original version was 9 September 2016 
#          revamped version included into PyFluxPro 
#          (see https://github.com/OzFlux/PyFluxPro for details)
#          new standalone version 26 September 2018
#                                                     Version: 1.0
#
#                                              2018, Cacilia Ewenz
#
# GUI to run Kljun etal 2015 or Korman and Meixner 2001 footprint
# or create windrose climatologies using OzFlux netcdf file format.
# INPUT: KM01  ==> e.g. footprint.txt
#        Kljun ==> e.g. footprint.txt
#
# CHANGES:
#
# 2018/09/21 - CME: adjusted to pfp_* scripts, Peter's io version and 
#                   date/time for Day/Month/Year climatology
# 2018/09/24 - CME: streamlined scripts from PyFluxPro to only what is necessary.
#                   original PyFluxPro scripts may be modified to include only 
#                   whats used in footprint calculations, original files are:
#    cfg.py, constants.py, meteorologicalfuntions.py pfp_io.py, pfp_log.py, pfp_utils.py
# ===================================================================
import ast
import copy
import datetime
#import logging
import matplotlib
matplotlib.use('TkAgg')
import numpy
import ntpath
import time
import Tkinter as tk
import tkMessageBox
from PIL import Image, ImageTk
import os
import sys
# The Lindsay Trap: check the scripts directory is present
if not os.path.exists("./scripts/"):
    print "PyFluxPro: the scripts directory is missing"
    sys.exit()
# since the scripts directory is there, try importing the modules
sys.path.append('./scripts')
import footprint_io
import footprint_log
# load the footprint routines ============
import footprint_fp
# load the windrose routines ============
import footprint_wr
# ========================================
# now check the logfiles and plots directories are present
dir_list = ["./logfiles/", "./plots/"]
for item in dir_list:
    if not os.path.exists(item):
        os.makedirs(item)
# start a log file with the current date and time in the name
t = time.localtime()
rundatetime = datetime.datetime(t[0], t[1], t[2], t[3], t[4], t[5]).strftime("%Y%m%d%H%M")
log_filename = 'footprint_'+rundatetime+'.log'
logger = footprint_log.init_logger(logger_name="footprint_log", file_handler=log_filename)


class qcgui(tk.Tk):
    """
        """
    def __init__(self, parent):
        tk.Tk.__init__(self,parent)
        self.parent = parent
        self.initialise()

    def option_not_implemented(self):
        self.do_progress(text='Option not implemented yet ...')
        logger.info(' Option not implemented yet ...')
    
    def initialise(self):
        self.org_frame = tk.Frame(self,height=600,width=500)
        self.org_frame.grid_propagate(0)

        self.org_frame.grid()

        im0 = Image.open('logo.png')
        tkimage0 = ImageTk.PhotoImage(im0)
        myvar0=tk.Label(self.org_frame,image = tkimage0)
        myvar0.image = tkimage0 # keep a reference
        myvar0.place(x=0, y=0, relwidth=1, relheight=1)

        #im1 = Image.open('ozflux_logo.gif')
        #tkimage1 = ImageTk.PhotoImage(im1)
        #myvar1=tk.Label(self.org_frame,image = tkimage1)
        #myvar1.image = tkimage1 # keep a reference
        #myvar1.place(x=280, y=0, width=220, height=63)

        #im2 = Image.open('tern_logo.jpg')
        #tkimage2 = ImageTk.PhotoImage(im2)
        #myvar2=tk.Label(self.org_frame,image = tkimage2)
        #myvar2.image = tkimage2 # keep a reference
        #myvar2.place(x=280, y=63, width=220, height=92)

        L1Label = tk.Label(self.org_frame,text='Kormann and Meixner, 2001')
        L1Label.place(x=15,y=525,width=190, height=20)
        doL1Button = tk.Button (self.org_frame, text="KM01", command=lambda:self.do_footprint(mode="kormei") )
        doL1Button.place(x=250,y=525,width=60, height=20)
        L2Label = tk.Label(self.org_frame,text='Kljun et al., 2015')
        L2Label.place(x=35,y=550,width=150, height=20)
        doL2Button = tk.Button (self.org_frame, text="Kljun", command=lambda:self.do_footprint(mode="kljun") )
        doL2Button.place(x=250,y=550,width=60, height=20)
        L3Label = tk.Label(self.org_frame,text='Windrose Climatology')
        L3Label.place(x=30,y=575,width=160, height=20)
        doL3Button = tk.Button (self.org_frame, text="WindRose", command=self.do_plotwindroseclimatology )
        doL3Button.place(x=250,y=575,width=70, height=20)

        closeplotwindowsButton = tk.Button (self.org_frame, text="Close plot windows", command=self.do_closeplotwindows )
        closeplotwindowsButton.place(x=370,y=550,width=120, height=20)

        quitButton = tk.Button (self.org_frame, text='Quit', command=self.do_quit )
        quitButton.place(x=440,y=575,width=50, height=20)
        self.progress = tk.Label(self.org_frame, text='Waiting for input ...')
        self.progress.grid(row=0,column=0,columnspan=6,sticky="W")

        # now we put together the menu, "Help" menu
        menubar = tk.Menu(self)
        helpmenu = tk.Menu(menubar,tearoff=0)
        helpmenu.add_command(label="Contents",command=self.do_helpcontents)
        helpmenu.add_command(label="About",command=self.option_not_implemented)
        menubar.add_cascade(label="Help",menu=helpmenu)
        self.config(menu=menubar)

    def do_closeplotwindows(self):
        """
            Close plot windows
            """
        self.do_progress(text='Closing plot windows ...')             # tell the user what we're doing
        logger.info(' Closing plot windows ...')
        matplotlib.pyplot.close('all')
        #fig_numbers = [n.num for n in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
        ##logger.info('  Closing plot windows: '+str(fig_numbers))
        #for n in fig_numbers:
            #matplotlib.pyplot.close(n)
        self.do_progress(text='Waiting for input ...')             # tell the user what we're doing
        logger.info(' Waiting for input ...')

    # ### footprint start
    def do_footprint(self,mode="kljun"):
        """
        Calls footprint_fp.footprint_main
        kljun  = Calculates the Kljun et al., 2015 footprint climatology.
        kormei = Calculates the Korman&Meixner, 2001 footprint climatology.
        """
        logger.info(' Starting footprint climatology')
        self.do_progress(text=' Calculating footprint climatology')
        if mode=="kljun":
            self.do_progress(text='Loading Kljun et al. control file ...')
            cf = footprint_io.load_controlfile(path='controlfiles')
            if len(cf)==0:
                self.do_progress(text='Waiting for input ...')
                return
        elif mode=="kormei":
            self.do_progress(text='Loading Kormann & Meixner control file ...')
            cf = footprint_io.load_controlfile(path='controlfiles')
            if len(cf)==0:
                self.do_progress(text='Waiting for input ...')
                return
        self.do_progress(text='Doing the '+mode+' footprint climatology')
        if "Options" not in cf:
            cf["Options"]={}
        cf["Options"]["call_mode"] = "interactive"
        footprint_fp.footprint_main(cf,mode)
        self.do_progress(text='Finished calculating footprint')
        logger.info(' Finished calculating footprint')
        logger.info("")
    # ### footprint end

    # === plot windroses
    def do_plotwindroseclimatology(self):
        """
           Plot windrose climatology. Set annual, monthly, or special timestep
        """
        cf = footprint_io.load_controlfile(path='controlfiles')
        if len(cf)==0:
            self.do_progress(text='Waiting for input ...')
            return
        if "Options" not in cf: cf["Options"]={}
        cf["Options"]["call_mode"] = "interactive"
        wrfilename = footprint_io.get_outfilenamefromcf(cf)
        footprint_wr.windrose_main(cf)
        self.do_progress(text='Finished calculating windrose climatology')
        logger.info(' Finished calculating windrose climatology')
        logger.info("")
        # ===

    def do_helpcontents(self):
        tkMessageBox.showinfo("Obi Wan says ...","Read the source, Luke!")

    def do_progress(self,text):
        """
            Update progress message in QC Data GUI
            """
        self.progress.destroy()
        self.progress = tk.Label(self.org_frame, text=text)
        self.progress.grid(row=11,column=0,columnspan=6,sticky="W")
        self.update()

    def do_closeplotwindows(self):
        """
            Close plot windows
            """
        self.do_progress(text='Closing plot windows ...')             # tell the user what we're doing
        logger.info(' Closing plot windows ...')
        matplotlib.pyplot.close('all')
        #fig_numbers = [n.num for n in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
        ##logger.info('  Closing plot windows: '+str(fig_numbers))
        #for n in fig_numbers:
            #matplotlib.pyplot.close(n)
        self.do_progress(text='Waiting for input ...')             # tell the user what we're doing
        logger.info(' Waiting for input ...')

    def do_quit(self):
        """
            Close plot windows and quit QC Data GUI
            """
        self.do_progress(text='Closing plot windows ...')             # tell the user what we're doing
        logger.info(' Closing plot windows ...')
        matplotlib.pyplot.close('all')
        self.do_progress(text='Quitting ...')                         # tell the user what we're doing
        logger.info(' Quitting ...')
        self.quit()

    def update_startenddate(self,startstr,endstr):
        """
            Read start and end timestamps from data and report in QC Data GUI
            """
        self.filestartValue.destroy()
        self.fileendValue.destroy()
        self.filestartValue = tk.Label(self.org_frame,text=startstr)
        self.filestartValue.grid(row=3,column=3,columnspan=3)
        self.fileendValue = tk.Label(self.org_frame,text=endstr)
        self.fileendValue.grid(row=4,column=3,columnspan=3)
        self.update()
    
if __name__ == "__main__":
    #log = qcutils.startlog('qc','logfiles/qc.log')
    qcGUI = qcgui(None)
    #main_title = cfg.version_name+' Main GUI '+cfg.version_number
    main_title = "TERN OzFlux FootPrint calulation V1.0"
    qcGUI.title(main_title)
    qcGUI.mainloop()
    qcGUI.destroy()

    logger.info('QC: All done')

