#!/usr/bin/env python

# Copyright 2012 by Hartmut Stoecker
# Contact: hartmut.stoecker@physik.tu-freiberg.de
#
# MergeConvert provides a graphical user interface to display X-ray diffraction data and some special tools for X-ray reflectivity measurements.
#
# The present version is 0.5.

import os, wx, sys
from wx.lib.mixins.listctrl import CheckListCtrlMixin

import matplotlib
matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.widgets import RectangleSelector
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas, NavigationToolbar2WxAgg as NavigationToolbar
    
from pylab import arange, argwhere, array, exp, float_, int32, float32, float64, fromfile, loadtxt, log, mod, pi, savetxt, sin, vstack, zeros
from StringIO import StringIO
from time import gmtime, strftime
from scipy.optimize import leastsq


class CheckListCtrl(wx.ListCtrl, CheckListCtrlMixin):
    def __init__(self, parent):
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT)
        CheckListCtrlMixin.__init__(self)

        
class MergeWindow(wx.Dialog):
    """ Dialog to select merge parameters. """

    def __init__(self, parent, id, title):
        wx.Dialog.__init__(self, parent, id, title)
        
        self.value = 0.0
        self.angle = array([0, 0])
        self.int = array([1, 1])
        self.a = app.frame.list.GetFirstSelected()
        self.b = app.frame.list.GetNextSelected(self.a)

        self.panel = wx.Panel(self)
        
        self.faktorlabel = wx.StaticText(self.panel, -1, "Faktor:")
        self.faktor = wx.TextCtrl(self.panel, -1, '0.0', style=wx.TE_PROCESS_ENTER)    
        self.faktor.Bind(wx.EVT_TEXT_ENTER, self.do_merge)
        self.faktor.Bind(wx.EVT_TEXT, self.on_action)
        
        self.cb_auto = wx.CheckBox(self.panel, -1, "Auto", style=wx.ALIGN_LEFT)
        self.cb_auto.SetValue(1)
        self.Bind(wx.EVT_CHECKBOX, self.do_merge, self.cb_auto)
                
        self.redrawbutton = wx.Button(self.panel, -1, "Neu zeichnen")
        self.Bind(wx.EVT_BUTTON, self.do_merge, self.redrawbutton)
        
        nametext = app.frame.oldname[self.a] + ' + ' + app.frame.oldname[self.b]
        self.namelabel = wx.StaticText(self.panel, -1, "Name:")
        self.name = wx.TextCtrl(self.panel, -1, nametext, style=wx.TE_PROCESS_ENTER, size=(240, 20))  
        self.name.Bind(wx.EVT_TEXT_ENTER, self.do_merge)
        
        self.okButton = wx.Button(self.panel, wx.ID_OK, 'OK', size=(90, 25))
        self.Bind(wx.EVT_BUTTON, self.on_OK, self.okButton)
        self.SetAffirmativeId(wx.ID_OK)

        self.closeButton = wx.Button(self.panel, wx.ID_CANCEL, 'Abbrechen', size=(90, 25))
        self.SetEscapeId(wx.ID_CANCEL)
        
        flags = wx.ALIGN_CENTER | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox1.Add(self.namelabel, 0, border=3, flag=flags)
        self.hbox1.Add(self.name, 0, border=3, flag=flags)
        
        self.hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox2.Add(self.faktorlabel, 0, border=3, flag=flags)
        self.hbox2.Add(self.faktor, 0, border=3, flag=flags)
        self.hbox2.Add(self.cb_auto, 0, border=3, flag=flags)
        self.hbox2.Add(self.redrawbutton, 0, border=3, flag=flags)
              
        self.hbox3 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox3.Add(self.okButton, 0, border=3, flag=flags)
        self.hbox3.Add(self.closeButton, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox1, 0, flag=wx.ALIGN_CENTER | wx.TOP | wx.EXPAND)
        self.vbox.Add(self.hbox2, 0, flag=wx.ALIGN_CENTER | wx.TOP | wx.EXPAND)
        self.vbox.AddSpacer(10)
        self.vbox.Add(self.hbox3, 0, flag=wx.ALIGN_CENTER | wx.TOP)
        
        self.panel.SetSizer(self.vbox)
        self.vbox.Fit(self)
        
        self.do_merge('', 1)
        
    def on_action(self, event):
        raw_value = self.faktor.GetValue().strip()
        try:
            self.value = float(raw_value)
        except ValueError:
            self.faktor.ChangeValue(str(self.value))
            
    def on_OK(self, event):
        self.do_merge(event)
        self.SetReturnCode(wx.ID_OK)
        self.Destroy()
            
    def do_merge(self, event, append=0):
        angle1 = app.frame.data[self.a][:,0]
        angle2 = app.frame.data[self.b][:,0]
        int1 = app.frame.data[self.a][:,1]
        int2 = app.frame.data[self.b][:,1]
        success = 0

        if angle1[0] < angle2[0] and angle1[-1] > angle2[0]:
        
            start1 = max( argwhere(angle1 >= angle2[0])[0,0], argwhere(int2 <= 1e5)[0,0] )
            start2 = argwhere(angle2.round(4) >= angle1[start1].round(4))[0,0]
            end1 = min( len(angle1), start1+len(angle2) )
            end2 = start2 - start1 + end1
            
            if self.cb_auto.IsChecked():
                fact = leastsq(lambda p: int1[start1:end1] - p*int2[start2:end2], 0.01)[0][0]
                self.value = fact
                self.faktor.ChangeValue(str(fact))
            else:
                fact = self.value
            
            diff = start1 - start2
            length = len(angle2) + diff
            angle = zeros(length)
            int = zeros(length)
            
            for i in arange(length):
                if i < start1:
                    angle[i] = angle1[i]
                    int[i] = int1[i]
                elif i <= start1+20 and i < end1:
                    angle[i] = angle1[i]
                    int[i] = (int1[i] + fact * int2[i-diff]) / 2
                else:
                    angle[i] = angle2[i-diff]
                    int[i] = fact * int2[i-diff]
                    
            success = 1
                    
        elif angle2[0] < angle1[0] and angle2[-1] > angle1[0]:
        
            start2 = max( argwhere(angle2 >= angle1[0])[0,0], argwhere(int1 <= 1e5)[0,0] )
            start1 = argwhere(angle1.round(4) >= angle2[start2].round(4))[0,0]
            end2 = min( len(angle2), start2+len(angle1) )
            end1 = start1 - start2 + end2
            
            if self.cb_auto.IsChecked():
                fact = leastsq(lambda p: int1[start1:end1] - p*int2[start2:end2], 0.01)[0][0]
                self.value = fact
                self.faktor.ChangeValue(str(fact))
            else:
                fact = self.value
            
            diff = start2 - start1
            length = len(angle1) + diff
            angle = zeros(length)
            int = zeros(length)
            
            for i in arange(length):
                if i < start2:
                    angle[i] = angle2[i]
                    int[i] = fact * int2[i]
                elif i <= start2+20 and i < end2:
                    angle[i] = angle2[i]
                    int[i] = (int1[i-diff] + fact * int2[i]) / 2
                else:
                    angle[i] = angle1[i-diff]
                    int[i] = int1[i-diff]
                    
            success = 1
                    
        else:
            fact = 0.1
            self.value = fact
            self.faktor.ChangeValue(str(fact))        
            angle = array([0, 0])
            int = array([1, 1])
            
            self.SetReturnCode(wx.ID_CANCEL)
            self.Destroy()
            wx.MessageBox('Verbinden dieser Daten nicht erfolgreich!','Fehler')
        
        if success:
            self.angle = angle
            self.int = int
            
            if append:
                self.append_data()
            else:
                self.replace_data()
            
            app.frame.update_list()
            app.frame.draw_figure()
            
    def append_data(self):
        app.frame.colors.append(app.frame.defaultcolors[mod(len(app.frame.newname), len(app.frame.defaultcolors))])
        app.frame.newname.append(os.path.splitext(self.name.GetValue())[0])
        app.frame.oldname.append(self.name.GetValue())
        app.frame.comment.append(app.frame.oldname[self.a] + ' + ' + app.frame.oldname[self.b] + ' (skaliert auf ' + str(self.value) + ')')
        app.frame.date.append(strftime("%d-%b-%Y, %H:%M:%S"))
        app.frame.wavelength.append(app.frame.wavelength[self.a])
        app.frame.omega.append(app.frame.omega[self.a])
        app.frame.twotheta.append(app.frame.twotheta[self.a])
        app.frame.scantype.append(app.frame.scantype[self.a])
        app.frame.scanaxis.append(app.frame.scanaxis[self.a])
        app.frame.first.append(self.angle[0])
        app.frame.range.append(self.angle[-1] - self.angle[0])
        app.frame.step.append(self.angle[1] - self.angle[0])
        app.frame.time.append(min(app.frame.time[self.a], app.frame.time[self.b]))
        app.frame.number.append(len(self.angle))
        app.frame.data.append(vstack((self.angle, self.int)).T)
        app.frame.checked.append(1)
        
    def replace_data(self):
        i = len(app.frame.newname) - 1
    
        app.frame.newname[i] = os.path.splitext(self.name.GetValue())[0]
        app.frame.oldname[i] = self.name.GetValue()
        app.frame.comment[i] = app.frame.oldname[self.a] + ' + ' + app.frame.oldname[self.b] + ' (skaliert auf ' + str(self.value) + ')'
        app.frame.date[i] = strftime("%d-%b-%Y, %H:%M:%S")
        app.frame.first[i] = self.angle[0]
        app.frame.range[i] = self.angle[-1] - self.angle[0]
        app.frame.step[i] = self.angle[1] - self.angle[0]
        app.frame.number[i] = len(self.angle)
        app.frame.data[i] = vstack((self.angle, self.int)).T
        app.frame.checked[i] = 1
        

class CorrectWindow(wx.Dialog):
    """ Dialog to select XRR correction parameters. """

    def __init__(self, parent, id, title):
        wx.Dialog.__init__(self, parent, id, title)
        
        self.beamsize = 0.1
        self.samplesize = 5.0
        self.a = app.frame.list.GetFirstSelected()
        self.angle = app.frame.data[self.a][:,0]
        self.int = app.frame.data[self.a][:,1]
        self.res = 1.0 * self.int

        self.panel = wx.Panel(self)
        
        self.bslabel = wx.StaticText(self.panel, -1, "Strahl (mm):")
        self.bs = wx.TextCtrl(self.panel, -1, str(self.beamsize), style=wx.TE_PROCESS_ENTER, size=(80, 20))    
        self.bs.Bind(wx.EVT_TEXT_ENTER, self.do_correct)
        self.bs.Bind(wx.EVT_TEXT, self.on_bs)
        
        self.sslabel = wx.StaticText(self.panel, -1, "Probe (mm):")
        self.ss = wx.TextCtrl(self.panel, -1, str(self.samplesize), style=wx.TE_PROCESS_ENTER, size=(80, 20))    
        self.ss.Bind(wx.EVT_TEXT_ENTER, self.do_correct)
        self.ss.Bind(wx.EVT_TEXT, self.on_ss)
        
        self.cb_tt = wx.CheckBox(self.panel, -1, "Winkel ist 2theta", style=wx.ALIGN_LEFT)
        self.cb_tt.SetValue(1)
        self.Bind(wx.EVT_CHECKBOX, self.do_correct, self.cb_tt)
        
        self.redrawbutton = wx.Button(self.panel, -1, "Neu zeichnen")
        self.Bind(wx.EVT_BUTTON, self.do_correct, self.redrawbutton)
        
        nametext = os.path.splitext(app.frame.oldname[self.a])[0] + '_corr'
        self.namelabel = wx.StaticText(self.panel, -1, "Name:")
        self.name = wx.TextCtrl(self.panel, -1, nametext, style=wx.TE_PROCESS_ENTER, size=(260, 20))  
        self.name.Bind(wx.EVT_TEXT_ENTER, self.do_correct)
        
        self.okButton = wx.Button(self.panel, wx.ID_OK, 'OK', size=(90, 25))
        self.Bind(wx.EVT_BUTTON, self.on_OK, self.okButton)
        self.SetAffirmativeId(wx.ID_OK)

        self.closeButton = wx.Button(self.panel, wx.ID_CANCEL, 'Abbrechen', size=(90, 25))
        self.SetEscapeId(wx.ID_CANCEL)
        
        flags = wx.ALIGN_CENTER | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox1.Add(self.namelabel, 0, border=3, flag=flags)
        self.hbox1.Add(self.name, 0, border=3, flag=flags)
        
        self.hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox2.Add(self.bslabel, 0, border=3, flag=flags)
        self.hbox2.Add(self.bs, 0, border=3, flag=flags)
        self.hbox2.Add(self.sslabel, 0, border=3, flag=flags)
        self.hbox2.Add(self.ss, 0, border=3, flag=flags)
        
        self.hbox3 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox3.Add(self.cb_tt, 0, border=3, flag=flags)
        self.hbox3.Add((110, 10), 0)
        self.hbox3.Add(self.redrawbutton, 0, border=3, flag=flags)
              
        self.hbox4 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox4.Add(self.okButton, 0, border=3, flag=flags)
        self.hbox4.Add(self.closeButton, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox1, 0, flag=wx.ALIGN_CENTER | wx.TOP | wx.EXPAND)
        self.vbox.Add(self.hbox2, 0, flag=wx.ALIGN_CENTER | wx.TOP | wx.EXPAND)
        self.vbox.Add(self.hbox3, 0, flag=wx.ALIGN_CENTER | wx.TOP | wx.EXPAND)
        self.vbox.AddSpacer(10)
        self.vbox.Add(self.hbox4, 0, flag=wx.ALIGN_CENTER | wx.TOP)
        
        self.panel.SetSizer(self.vbox)
        self.vbox.Fit(self)
        
        self.do_correct('', 1)
        
    def on_bs(self, event):
        raw_value = self.bs.GetValue().strip()
        try:
            self.beamsize = float(raw_value)
        except ValueError:
            self.bs.ChangeValue(str(self.beamsize))
            
    def on_ss(self, event):
        raw_value = self.ss.GetValue().strip()
        try:
            self.samplesize = float(raw_value)
        except ValueError:
            self.ss.ChangeValue(str(self.samplesize))
            
    def on_OK(self, event):
        self.do_correct(event)
        self.SetReturnCode(wx.ID_OK)
        self.Destroy()
            
    def do_correct(self, event, append=0):
        if self.cb_tt.IsChecked():
            omega = self.angle / 2.0
        else:
            omega = self.angle
    
        fakt = self.samplesize / self.beamsize * abs( sin( omega*pi/180 ) )
        fakt = fakt * (fakt<1) + 1.0 * (fakt>1)
        self.res = self.int / fakt
                 
        if append:
            self.append_data()
        else:
            self.replace_data()
        
        app.frame.update_list()
        app.frame.draw_figure()
            
    def append_data(self):
        app.frame.colors.append(app.frame.defaultcolors[mod(len(app.frame.newname), len(app.frame.defaultcolors))])
        app.frame.newname.append(os.path.splitext(self.name.GetValue())[0])
        app.frame.oldname.append(self.name.GetValue())
        app.frame.comment.append(app.frame.oldname[self.a] + ' (korrigiert mit Strahl = ' + str(self.beamsize) + ' mm und Probe = ' + str(self.samplesize) + ' mm)')
        app.frame.date.append(strftime("%d-%b-%Y, %H:%M:%S"))
        app.frame.wavelength.append(app.frame.wavelength[self.a])
        app.frame.omega.append(app.frame.omega[self.a])
        app.frame.twotheta.append(app.frame.twotheta[self.a])
        app.frame.scantype.append(app.frame.scantype[self.a])
        app.frame.scanaxis.append(app.frame.scanaxis[self.a])
        app.frame.first.append(self.angle[0])
        app.frame.range.append(self.angle[-1] - self.angle[0])
        app.frame.step.append(self.angle[1] - self.angle[0])
        app.frame.time.append(app.frame.time[self.a])
        app.frame.number.append(len(self.angle))
        app.frame.data.append(vstack((self.angle, self.res)).T)
        app.frame.checked.append(1)
        
    def replace_data(self):
        i = len(app.frame.newname) - 1
    
        app.frame.newname[i] = os.path.splitext(self.name.GetValue())[0]
        app.frame.oldname[i] = self.name.GetValue()
        app.frame.comment[i] = app.frame.oldname[self.a] + ' (korrigiert mit Strahl = ' + str(self.beamsize) + ' mm und Probe = ' + str(self.samplesize) + ' mm)'
        app.frame.date[i] = strftime("%d-%b-%Y, %H:%M:%S")
        app.frame.first[i] = self.angle[0]
        app.frame.range[i] = self.angle[-1] - self.angle[0]
        app.frame.step[i] = self.angle[1] - self.angle[0]
        app.frame.number[i] = len(self.angle)
        app.frame.data[i] = vstack((self.angle, self.res)).T
        app.frame.checked[i] = 1
        

class MainFrame(wx.Frame):
    """ The main frame of the application."""
    
    def __init__(self):
        wx.Frame.__init__(self, None, -1, 'MergeConvert - Anzeige, Konvertierung und Verbinden von Diffraktometer-Dateien')
        
        self.newname = []
        self.oldname = []
        self.comment = []
        self.date = []
        self.wavelength = []
        self.omega = []
        self.twotheta = []
        self.scantype = []
        self.scanaxis = []
        self.first = []
        self.range = []
        self.step = []
        self.time = []
        self.number = []
        self.data = []
        
        self.savename = ''
        self.checked = []
        self.redraw = 1
        self.defaultcolors = [(0,0,1), (0,0.5,0), (1,0,0), (0,0.75,0.75), (0.75,0,0.75), (0.75,0.75,0), (0,0,0)]
        self.colors = []
        
        self.create_menu()
        self.create_status_bar()
        self.create_main_panel()
        self.draw_figure()

    def create_menu(self):
        self.menubar = wx.MenuBar()
        
        menu_file = wx.Menu()
        m_load = menu_file.Append(-1, "&Datei laden...\tCtrl-O", "Daten aus Datei laden")
        self.Bind(wx.EVT_MENU, self.on_open_file, m_load)
        m_delete = menu_file.Append(-1, "&Datei entfernen...\tCtrl-L", "Bestehende Datei entfernen")
        self.Bind(wx.EVT_MENU, self.on_delete_file, m_delete)
        m_color = menu_file.Append(-1, "&Farbe wechseln...\tCtrl-F", "Farbe wechseln")
        self.Bind(wx.EVT_MENU, self.on_color, m_color)
        m_merge = menu_file.Append(-1, "&Daten verbinden...\tCtrl-M", "Daten verbinden")
        self.Bind(wx.EVT_MENU, self.on_merge, m_merge)
        m_correct = menu_file.Append(-1, "&XRR-Korrektur...\tCtrl-K", "XRR-Daten korrigieren")
        self.Bind(wx.EVT_MENU, self.on_correct, m_correct)
        menu_file.AppendSeparator()
        m_savetext = menu_file.Append(-1, "&Daten speichern...\tCtrl-S", "Daten als Textdatei speichern")
        self.Bind(wx.EVT_MENU, self.on_save_text, m_savetext)
        m_saveplot = menu_file.Append(-1, "&Grafik speichern...\tCtrl-G", "Grafik als Datei speichern")
        self.Bind(wx.EVT_MENU, self.on_save_plot, m_saveplot)
        menu_file.AppendSeparator()
        m_exit = menu_file.Append(-1, "&Beenden\tCtrl-X", "Programm verlassen")
        self.Bind(wx.EVT_MENU, self.on_exit, m_exit)
        
        menu_help = wx.Menu()
        m_about = menu_help.Append(-1, "&Hilfe\tF1", "Hilfe zum Programm")
        self.Bind(wx.EVT_MENU, self.on_about, m_about)
        
        self.menubar.Append(menu_file, "&Datei")
        self.menubar.Append(menu_help, "&Hilfe")
        self.SetMenuBar(self.menubar)

    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas 
             * mpl navigation toolbar
             * Control panel for interaction
        """
        
        self.spWindow = wx.SplitterWindow(self)
        self.panel = wx.Panel(self.spWindow)
        self.graphpanel = wx.Panel(self.spWindow)
        self.spWindow.SetMinimumPaneSize(100)
        self.spWindow.SplitHorizontally(self.panel, self.graphpanel, 150)
        
        # Create the mpl Figure and FigCanvas objects: 9x6 inches, 100 dots-per-inch
        self.dpi = 100
        self.fig = Figure((9.1, 6), dpi=self.dpi)
        self.canvas = FigCanvas(self.graphpanel, -1, self.fig)
        self.canvas.mpl_connect('motion_notify_event', self.on_UpdateCursor)
        self.canvas.mpl_connect('resize_event', self.on_Resize)
        self.canvas.mpl_connect('button_press_event', self.on_Press)
        self.canvas.mpl_connect('scroll_event', self.on_Scroll)
        
        self.log = wx.TextCtrl(self.graphpanel, -1, size=(666,66), style = wx.TE_MULTILINE | wx.TE_READONLY | wx.TE_RICH2)
        self.log.SetFont(wx.Font(8, wx.MODERN, wx.NORMAL, wx.NORMAL, 0, "Sans"))
        sys.stdout = RedirectText(self.log, 'black')
        sys.stderr = RedirectText(self.log, 'red')
        
        self.vbox_g = wx.BoxSizer(wx.VERTICAL)
        self.vbox_g.Add(self.canvas, 1, wx.EXPAND)        
        self.vbox_g.Add(self.log, 0, wx.EXPAND)
        
        self.graphpanel.SetSizer(self.vbox_g)
        self.vbox_g.Fit(self)
        
        # Create the top panel
        self.loadbutton = wx.Button(self.panel, -1, "Datei\nladen", size=(70,35))
        self.Bind(wx.EVT_BUTTON, self.on_open_file, self.loadbutton)
              
        self.deletebutton = wx.Button(self.panel, -1, "Datei\nentfernen", size=(70,35))
        self.Bind(wx.EVT_BUTTON, self.on_delete_file, self.deletebutton)
        
        self.colorbutton = wx.Button(self.panel, -1, "Farbe\nwechseln", size=(70,35))
        self.Bind(wx.EVT_BUTTON, self.on_color, self.colorbutton)
        
        self.mergebutton = wx.Button(self.panel, -1, "Daten\nverbinden", size=(70,35))
        self.Bind(wx.EVT_BUTTON, self.on_merge, self.mergebutton)
        
        self.correctbutton = wx.Button(self.panel, -1, "XRR-\nKorrektur", size=(70,35))
        self.Bind(wx.EVT_BUTTON, self.on_correct, self.correctbutton)

        self.savetextbutton = wx.Button(self.panel, -1, "Daten\nspeichern", size=(70,35))
        self.Bind(wx.EVT_BUTTON, self.on_save_text, self.savetextbutton)
        
        self.saveplotbutton = wx.Button(self.panel, -1, "Grafik\nspeichern", size=(70,35))
        self.Bind(wx.EVT_BUTTON, self.on_save_plot, self.saveplotbutton)
        
        self.cb_cps = wx.CheckBox(self.panel, -1, "Counts pro Sekunde", style=wx.ALIGN_LEFT)
        self.cb_cps.SetValue(1)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb, self.cb_cps)
        
        self.cb_log = wx.CheckBox(self.panel, -1, "y logarithmisch", style=wx.ALIGN_LEFT)
        self.cb_log.SetValue(1)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb, self.cb_log)
        
        self.cb_grid = wx.CheckBox(self.panel, -1, "Gitternetz anzeigen", style=wx.ALIGN_LEFT)
        self.cb_grid.SetValue(1)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb, self.cb_grid)

        self.cb_header = wx.CheckBox(self.panel, -1, "Kopfdaten mitspeichern", style=wx.ALIGN_LEFT)
        self.cb_header.SetValue(1)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb, self.cb_header)
        
        self.cb_legend = wx.CheckBox(self.panel, -1, "Legende anzeigen", style=wx.ALIGN_LEFT)
        self.cb_legend.SetValue(1)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb, self.cb_legend)
        
        self.cb_showlog = wx.CheckBox(self.panel, -1, "Log anzeigen", style=wx.ALIGN_LEFT)
        self.cb_showlog.SetValue(1)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_showlog, self.cb_showlog)
        
        self.list = CheckListCtrl(self.panel)
        self.list.OnCheckItem = self.on_CheckItem
        self.list.InsertColumn(0, "Dateiname", width=240)
        self.list.InsertColumn(1, "Startwinkel", width=70)
        self.list.InsertColumn(2, "Endwinkel", width=70)
        self.list.InsertColumn(3, "Schrittweite", width=75)
        self.list.InsertColumn(4, "Messzeit", width=65)
        self.list.InsertColumn(5, "Datum", width=120)
        self.list.InsertColumn(6, "Kommentar", width=240)
        
        # Layout with box sizers
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        
        self.vbox1 = wx.BoxSizer(wx.VERTICAL)
        self.vbox1.Add(self.cb_cps, 0, border=3, flag=flags)
        self.vbox1.Add(self.cb_log, 0, border=3, flag=flags)
        
        self.vbox2 = wx.BoxSizer(wx.VERTICAL)
        self.vbox2.Add(self.cb_grid, 0, border=3, flag=flags)
        self.vbox2.Add(self.cb_header, 0, border=3, flag=flags)
        
        self.vbox3 = wx.BoxSizer(wx.VERTICAL)
        self.vbox3.Add(self.cb_legend, 0, border=3, flag=flags)
        self.vbox3.Add(self.cb_showlog, 0, border=3, flag=flags)
 
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox.Add(self.loadbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.deletebutton, 0, border=3, flag=flags)
        self.hbox.Add(self.colorbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.mergebutton, 0, border=3, flag=flags)
        self.hbox.Add(self.correctbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.savetextbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.saveplotbutton, 0, border=3, flag=flags)
        self.hbox.AddSpacer(7)      
        self.hbox.Add(self.vbox1, 0, flag = wx.CENTER)
        self.hbox.Add(self.vbox2, 0, flag = wx.CENTER)
        self.hbox.Add(self.vbox3, 0, flag = wx.CENTER)
        
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT)
        self.vbox.Add(self.list, 1, border=3, flag = wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND)   
        
        self.panel.SetSizer(self.vbox)
        self.vbox.Fit(self.panel)
    
    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()
        self.statusbar.SetFieldsCount(3)
        self.statusbar.SetStatusWidths([-1,100,110])

    def draw_figure(self):
        """ Redraws the figure. """
        self.fig.clear()
        
        # Since we have only one plot, we can use add_axes instead of add_subplot, but then the subplot configuration tool in the navigation toolbar wouldn't work.
        # self.axes = self.fig.add_axes([0, 0, 1, 1])
        self.axes = self.fig.add_subplot(111)
        self.axes.locator_params(axis = 'x', nbins = 10)
        self.axes.grid(self.cb_grid.IsChecked())
        
        self.rs1 = RectangleSelector(self.axes, self.on_Zoom, drawtype='box', button=1, minspanx=5, minspany=5, spancoords='pixels')
        self.rs2 = RectangleSelector(self.axes, self.on_Move, drawtype='box', button=3)   # only 'box' will work without problems here
        
        checked = 0
        for i in arange(len(self.oldname)):
            if self.list.IsChecked(i):
                checked += 1
                if self.cb_cps.IsChecked() and self.time[i] != 0:
                    self.axes.plot(self.data[i][:,0], self.data[i][:,1]/self.time[i], label=self.oldname[i], color=self.colors[i])
                else:
                    self.axes.plot(self.data[i][:,0], self.data[i][:,1], label=self.oldname[i], color=self.colors[i])
        
        if self.cb_log.IsChecked() and self.data != []:
            self.axes.set_yscale('log')
        
        if self.cb_legend.IsChecked() and checked > 0:
            prop = matplotlib.font_manager.FontProperties(size=8) 
            self.axes.legend(loc=0, prop=prop)
        
        x = len(self.oldname)
        if x == 1:
            titeltext = 'Plot of X-ray diffraction data (' + str(checked) + ' of ' + str(x) + ' dataset selected)'
        else:
            titeltext = 'Plot of X-ray diffraction data (' + str(checked) + ' of ' + str(x) + ' datasets selected)'
        self.titel = self.axes.set_title(titeltext, fontsize=10)
               
        self.axes.tick_params(axis='both', labelsize=8)
        self.axes.set_xlabel('Angle (deg)', fontsize=10)
        
        if self.cb_cps.IsChecked():
            self.axes.set_ylabel('Intensity (cps)', fontsize=10)
        else:
            self.axes.set_ylabel('Intensity (counts)', fontsize=10)
        
        self.on_Resize('')
        self.canvas.draw()
        
    def on_Zoom(self, eclick, erelease):
        'eclick and erelease are the press and release events'
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata
        if x2 > x1:
            self.axes.set_xlim(min(x1,x2), max(x1,x2))
            self.axes.set_ylim(min(y1,y2), max(y1,y2))
        else:
            self.axes.autoscale()
        self.canvas.draw()
                
    def on_Move(self, eclick, erelease):
        'eclick and erelease are the press and release events'
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata
        xl, xm = self.axes.get_xlim()
        yl, ym = self.axes.get_ylim()
        
        self.axes.set_xlim(xl+x1-x2, xm+x1-x2)
        
        if self.axes.get_yscale() == 'log':
            self.axes.set_ylim(exp(log(yl)+log(y1)-log(y2)), exp(log(ym)+log(y1)-log(y2)))
        else:
            self.axes.set_ylim(yl+y1-y2, ym+y1-y2)
            
        self.canvas.draw()
        
    def on_Scroll(self, event):
        x, y = event.xdata, event.ydata
        xl, xm = self.axes.get_xlim()
        yl, ym = self.axes.get_ylim()
        factor = 1 - event.step / 20.0
        
        if not 'shift' in str(event.key):
            xln = x - factor * (x - xl)
            xmn = x + factor * (xm - x)
            self.axes.set_xlim(xln, xmn)
        if not 'control' in str(event.key):
            if self.axes.get_yscale() == 'log':
                yln_log = log(y) - factor * (log(y) - log(yl))
                ymn_log = log(y) + factor * (log(ym) - log(y))
                self.axes.set_ylim(exp(yln_log), exp(ymn_log))
            else:
                yln = y - factor * (y - yl)
                ymn = y + factor * (ym - y)
                self.axes.set_ylim(yln, ymn)
                
        self.canvas.draw()
            
    def on_Press(self, event):
        if hasattr(event, 'dblclick'):
            if event.dblclick:
                self.axes.autoscale()
                self.canvas.draw()
        else:
            print('Matplotlib Version 1.2 oder neuer zur Doppelklick-Erkennung erforderlich!')
        
    def on_UpdateCursor(self, event):
        if event.inaxes:
            if abs(event.xdata) > 1e-3 and abs(event.xdata) < 1e5:
                text1 = 'x = %.5f' %event.xdata
            else:   
                text1 = 'x = %.4e' %event.xdata
            if abs(event.ydata) > 1e-3 and abs(event.ydata) < 1e5:
                text2 = 'y = %.5f' %event.ydata
            else:   
                text2 = 'y = %.4e' %event.ydata
            self.statusbar.SetStatusText('Linksklick: Zoomen (links>rechts) bzw. Autozoom (rechts>links), Rechtsklick: Bewegen, Scrollen: Zoomen (Strg: x, Shift: y)', 0)
            self.statusbar.SetStatusText(text1, 1)
            self.statusbar.SetStatusText(text2, 2)
        else:
            self.statusbar.SetStatusText('', 0)
            self.statusbar.SetStatusText('', 1)
            self.statusbar.SetStatusText('', 2)
            
    def on_Resize(self, event):
        try:
            # self.fig.tight_layout(pad=1.0)
            x, y = self.fig.get_size_inches()
            self.fig.subplots_adjust(left=0.8/x, right=1-0.2/x, bottom=0.4/y, top=1-0.35/y)
            self.titel.set_position((0.5, 1+0.05/y))
        except Exception as error:
            return
            # print error
    
    def on_CheckItem(self, index, flag):
        if self.list.IsChecked(index):
            self.checked[index] = 1
        else:
            self.checked[index] = 0
        if self.redraw:
            self.draw_figure()
        
    def update_list(self):
        self.list.DeleteAllItems()
    
        for i in arange(len(self.oldname)):
            self.list.InsertStringItem(i, self.oldname[i])
            self.list.SetStringItem(i, 1, str(self.first[i]))
            self.list.SetStringItem(i, 2, str(self.range[i]+self.first[i]))
            self.list.SetStringItem(i, 3, str(self.step[i]))
            self.list.SetStringItem(i, 4, str(self.time[i]))
            self.list.SetStringItem(i, 5, self.date[i])
            self.list.SetStringItem(i, 6, self.comment[i])
        
        self.redraw = 0
        for i in arange(len(self.oldname)):
            if self.checked[i] == 1:
                self.list.CheckItem(i)
        self.redraw = 1
    
    def on_cb(self, event):
        self.draw_figure()
   
    def on_color(self, event):
        x = self.list.GetSelectedItemCount()
        
        if x == 1:
            i = self.list.GetFirstSelected()
            
            data = wx.ColourData()
            data.SetChooseFull(True)
            data.SetCustomColour(0, (10, 80, 161))
            data.SetCustomColour(1, (46, 144, 40))
            data.SetCustomColour(2, (178, 0, 38))
            data.SetCustomColour(3, (78, 188, 206))
            data.SetCustomColour(4, (230, 110, 1))        
            data.SetColour(wx.Colour(self.colors[i][0]*255, self.colors[i][1]*255, self.colors[i][2]*255))
            dlg = wx.ColourDialog(self, data)
           
            if dlg.ShowModal() == wx.ID_OK:
                res = dlg.GetColourData().Colour
                self.colors[i] = (res[0]/255.0, res[1]/255.0, res[2]/255.0)
                self.draw_figure()
                
        else:
            wx.MessageBox('Bitte genau einen Datensatz in der Liste markieren!','Fehler')
        
    def on_merge(self, event):
        x = self.list.GetSelectedItemCount()
        
        if x == 2:
            i = self.list.GetFirstSelected()
            j = self.list.GetNextSelected(i)
            
            if (self.step[i] - self.step[j]) < 1e-7:
                dlg = MergeWindow(None, -1, 'Daten verbinden')
                
                if dlg.ShowModal() == wx.ID_CANCEL:
                    i = len(self.newname) - 1
                    del self.colors[i]
                    del self.newname[i]
                    del self.oldname[i]
                    del self.comment[i]
                    del self.date[i]
                    del self.wavelength[i]
                    del self.omega[i]
                    del self.twotheta[i]
                    del self.scantype[i]
                    del self.scanaxis[i]
                    del self.first[i]
                    del self.range[i]
                    del self.step[i]
                    del self.time[i]
                    del self.number[i]
                    del self.data[i]
                    del self.checked[i]
                    self.update_list()
                    self.draw_figure()
                
            else:
                wx.MessageBox('Bitte nur Dateien mit gleicher Schrittweite verbinden!','Fehler')
            
        else:
            wx.MessageBox('Bitte genau zwei Dateien in der Liste markieren!','Fehler')
            
    def on_correct(self, event):
        x = self.list.GetSelectedItemCount()
        
        if x == 1:
            dlg = CorrectWindow(None, -1, 'XRR-Daten korrigieren')
                
            if dlg.ShowModal() == wx.ID_CANCEL:
                i = len(self.newname) - 1
                del self.colors[i]
                del self.newname[i]
                del self.oldname[i]
                del self.comment[i]
                del self.date[i]
                del self.wavelength[i]
                del self.omega[i]
                del self.twotheta[i]
                del self.scantype[i]
                del self.scanaxis[i]
                del self.first[i]
                del self.range[i]
                del self.step[i]
                del self.time[i]
                del self.number[i]
                del self.data[i]
                del self.checked[i]
                self.update_list()
                self.draw_figure()
            
        else:
            wx.MessageBox('Bitte genau einen Datensatz in der Liste markieren!','Fehler')
       
    def on_delete_file(self, event):
        x = self.list.GetSelectedItemCount()
        dodelete = 0
        
        if x == 0:
            wx.MessageBox('Bitte mindestens einen Datensatz in der Liste markieren!','Fehler')
        elif x == 1:
            dodelete = 1
        else:
            dlg = wx.MessageDialog(self, 'Alle ' + str(x) + ' markierten Dateien aus der Liste entfernen?', 'Frage', wx.YES_NO | wx.YES_DEFAULT | wx.ICON_QUESTION)
            if dlg.ShowModal() == wx.ID_YES:
                dodelete = 1
                
        if dodelete:
            for i in arange(len(self.checked)-1, -0.5, -1, dtype=int):
                if self.list.IsSelected(i):
                    del self.newname[i]
                    del self.oldname[i]
                    del self.comment[i]
                    del self.date[i]
                    del self.wavelength[i]
                    del self.omega[i]
                    del self.twotheta[i]
                    del self.scantype[i]
                    del self.scanaxis[i]
                    del self.first[i]
                    del self.range[i]
                    del self.step[i]
                    del self.time[i]
                    del self.number[i]
                    del self.data[i]
                    del self.checked[i]
                    del self.colors[i]
            self.update_list()
            self.draw_figure()
        
    def on_open_file(self, event):
        file_choices = "Daten-Typen (*.njc, *.raw, *.udf, *.x00, *.txt, *.dat)|*.njc;*.raw;*.udf;*.x00;*.txt;*.dat|Seifert (*.njc)|*.njc|Bruker (*.raw)|*.raw|Philips (*.udf)|*.udf|Philips (*.x00)|*.x00|TXT-Datei (*.txt)|*.txt|DAT-Datei (*.dat)|*.dat|Alle Dateien (*.*)|*.*"
        dlg = wx.FileDialog(self, "Datei laden", "", "", file_choices, wx.OPEN|wx.MULTIPLE)
        
        if dlg.ShowModal() == wx.ID_OK:
            filenames = dlg.GetFilenames()
            paths = dlg.GetPaths()
            
            for i in arange(len(paths)):
                filename = filenames[i]
                path = paths[i]
            
                oldname = filename
                newname, extension = os.path.splitext(filename)
                ext = extension.lstrip('.').lower()
                data = [[0,1],[0,1]]
                
                if ext == 'njc':
                    fd = open(path, 'rb')
                    z = fd.read(2000)
                    fd.close()
                    
                    comment1 = z[1008:1068].split('\x00')[0]
                    try:
                        offset1 = int(z[1324])
                    except:
                        offset1 = int(z[1328])
                    offset2 = ord(z[1340+offset1])
                    comment2_start = 1352 + offset1
                    comment2_end = comment2_start + offset2 - 9
                    comment2 = z[comment2_start:comment2_end]
                    
                    if comment1 != '' and comment2 != '':
                        comment = comment1 + ' - ' + comment2
                    else:
                        comment = comment1 + comment2
                    
                    offset = 1428 + offset1 + offset2
                    
                    fd = open(path, 'rb')
                    fd.read(offset)
                    data32 = fromfile(file=fd, dtype=float32)
                    fd.close()
                    fd = open(path, 'rb')
                    fd.read(offset)
                    data64 = fromfile(file=fd, dtype=float64)
                    fd.close()
                    
                    first = round(data64[20], 6)
                    range = round(data64[21], 6) - first
                    step = round(data32[44], 6)
                    number = int(round(range / step)) + 1
                    time = round(data32[46], 6)

                    angle = first + step * arange(number)
                    intens = data32[-40-number:-40] * time
                    data = vstack((angle, intens)).T
        
                    date = '???'
                    wavelength = str(data64[2])
                    omega = '???'
                    twotheta = '???'
                    scantype = '???'
                    scanaxis = '???'
                                   
                elif ext == 'raw':
                    fd = open(path, 'rb')
                    data_32 = fromfile(file=fd, dtype=float32)
                    fd.close()
                    fd = open(path, 'rb')
                    data_64 = fromfile(file=fd, dtype=float64)
                    fd.close()
                    fd = open(path, 'rb')
                    data_i = fromfile(file=fd, dtype=int32)
                    fd.close()
                    fd = open(path, 'rb')
                    data_s = fromfile(file=fd, dtype='a326')
                    fd.close()
                                                                     
                    comment = data_s[1][:160]
                    date = data_s[0][16:24] + ', ' + data_s[0][26:34]
                    wavelength = str(data_64[77])
                    omega = str(data_64[90])
                    twotheta = str(data_64[91])

                    if data_i[227] == 0:
                        scantype = 'Locked Coupled'
                        scanaxis = '2Theta'
                    elif data_i[227] == 1:
                        scantype = 'Unlocked Coupled'
                        scanaxis = '2Theta'
                    elif data_i[227] == 2:
                        scantype = 'Detector Scan'
                        scanaxis = '2Theta'
                    elif data_i[227] == 3:
                        scantype = 'Rocking Curve'
                        scanaxis = 'Omega'
                    elif data_i[227] == 9999:
                        scantype = 'Tube Scan'
                        scanaxis = 'Omega'
                    else:
                        scantype = '???'
                        scanaxis = '???'
                    
                    if scanaxis == 'Omega':
                        first = round(data_64[90], 6)
                    else:
                        first = round(data_64[91], 6)
                    
                    step = round(data_64[111], 6)
                    time = round(data_32[226], 6)
                    number = data_i[179]
                    range = step * (number - 1)
                    
                    angle = first + step * arange(number)
                    index = 254 + data_i[242] / 4
                    intens = data_32[index:index+number]          
                    data = vstack((angle, intens)).T
                    
                elif ext == 'udf':
                    f = open(path, 'r')
                    data_start = 0
                    intens = []
                    
                    for line in f:
                        if "SampleIdent" in line:
                            comment = line.split(',')[1]
                        if "FileDateTime" in line:
                            date = line.split(',')[1].lstrip(' ').replace(' ',', ')
                        if "LabdaAlpha1" in line:
                            l1 = float(line.split(',')[1])
                        if "LabdaAlpha2" in line:
                            l2 = float(line.split(',')[1])
                        if "RatioAlpha21" in line:
                            ratio = float(line.split(',')[1])
                        if "ScanType" in line:
                            scantype = line.split(',')[1].capitalize()
                        if "DataAngleRange" in line:
                            first = float(line.split(',')[1])
                            range = float(line.split(',')[2]) - first
                        if "ScanStepSize" in line:
                            step = float(line.split(',')[1])
                        if "ScanStepTime" in line:
                            time = float(line.split(',')[1])
                        if data_start == 1:
                            for x in line.rstrip('/\n').split(','):
                                if not x.strip() == '':
                                    intens.append(x.strip())
                        if "RawScan" in line:
                            data_start = 1
                    f.close()
                    
                    wavelength = str((l1 + l2*ratio) / (1 + ratio))
                    omega = '???'
                    twotheta = '???'
                    scanaxis = '???'
                    number = int(round(range / step)) + 1
                    
                    angle = first + step * arange(number)
                    data = vstack((angle, float_(intens))).T
                    
                elif ext == 'x00':
                    f = open(path, 'r')
                    header_length = 1
                    
                    for line in f:
                        # if "FileName" in line:
                            # oldname = line.split()[1]
                        if "Sample" in line:
                            comment = line.split()[1]
                        if "FileDateTime" in line:
                            date = line.split()[1].replace(',',', ')
                        if "Wavelength" in line:
                            wavelength = line.split()[1]
                        if "Omega," in line:
                            omega = line.split()[1]
                        if "TwoTheta" in line:
                            twotheta = line.split()[1]
                        if "ScanType" in line:
                            scantype = line.split()[1].capitalize()
                        if "ScanAxis" in line:
                            scanaxis = line.split()[1]
                        if "FirstAngle" in line:
                            first = float(line.split()[1])
                        if "ScanRange" in line:
                            range = float(line.split()[1])
                        if "StepWidth" in line:
                            step = float(line.split()[1])
                        if "TimePerStep" in line:
                            time = float(line.split()[1])
                        if "NrOfData" in line:
                            number = int(line.split()[1])
                        if "ScanData" in line:
                            break
                        else:
                            header_length += 1
                    f.close()
                    
                    angle = first + step * arange(number)
                    intens = loadtxt(path, skiprows=header_length) * time
                    data = vstack((angle, intens)).T
                
                else:
                    comment = '???'
                    date = '???'
                    wavelength = '???'
                    omega = '???'
                    twotheta = '???'
                    scantype = '???'
                    scanaxis = '???'
                    time = 0.0
                    first = 0.0
                    range = 0.0
                    step = 0.0
                    number = 0.0
                    
                    try:
                        use_cps = 0
                        header_length = 0
                        
                        f = open(path, 'r')
                        
                        for line in f:
                            splitting = line.lstrip(line.split(': ')[0] + ': ').strip()
                            # if "Original name:" in line:
                                # oldname = splitting
                            if "Comment:" in line:
                                comment = splitting
                            if "Original date:" in line:
                                date = splitting
                            if "Wavelength:" in line:
                                wavelength = splitting
                            if "Omega:" in line:
                                omega = splitting
                            if "2Theta:" in line:
                                twotheta = splitting
                            if "Scan type:" in line:
                                scantype = splitting
                            if "Scan axis:" in line:
                                scanaxis = splitting
                            if "Time per step:" in line:
                                time = float(splitting)
                            if "cps" in line:
                                use_cps = 1
                                
                            try:
                                data = loadtxt(path, skiprows=header_length)
                                break
                            except:
                                header_length += 1
                        
                        f.close()
                        
                        first = data[0,0]
                        range = data[-1,0] - data[0,0]
                        step = data[1,0] - data[0,0]
                        number = len(data[:,0])
                        
                        if use_cps and time != 0:
                            data[:,1] = data[:,1] * time
                            
                    except Exception as error:
                        wx.MessageBox('Beim Laden der Datei:\n\n' + path + '\n\ntrat folgender Fehler auf:\n\n' + str(error), 'Fehler beim Laden der Datei')
                
                if date == '???':
                    date = strftime("%d-%b-%Y, %H:%M:%S", gmtime(os.path.getmtime(path)))
                
                self.colors.append(self.defaultcolors[mod(len(self.newname), len(self.defaultcolors))])
                self.newname.append(newname)
                self.oldname.append(oldname)
                self.comment.append(comment)
                self.date.append(date)
                self.wavelength.append(wavelength)
                self.omega.append(omega)
                self.twotheta.append(twotheta)
                self.scantype.append(scantype)
                self.scanaxis.append(scanaxis)
                self.first.append(first)
                self.range.append(range)
                self.step.append(step)
                self.time.append(time)
                self.number.append(number)
                self.data.append(data)
                self.checked.append(1)
                            
                self.flash_status_message("%s geladen." % filename)
                self.update_list()
                self.draw_figure()

        dlg.Destroy()
        
    def on_save_plot(self, event):
        file_choices = "Portable Network Graphics (*.png)|*.png|Encapsulated Postscript (*.eps)|*.eps|Portable Document Format (*.pdf)|*.pdf|Scalable Vector Graphics (*.svg)|*.svg"
        filename = self.savename + ".png"
        dlg = wx.FileDialog(self, message="Grafik speichern unter...", defaultFile=filename, wildcard=file_choices, style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.savename = os.path.splitext(dlg.GetFilename())[0]
            self.canvas.print_figure(path, dpi=3*self.dpi)
            self.flash_status_message("Grafik gespeichert unter %s" % path)
        
    def on_save_text(self, event):
        i = self.list.GetFirstSelected()
        
        if i <= -1:
            wx.MessageBox('Bitte mindestens einen Datensatz in der Liste markieren!','Fehler')
        
        while i > -1:
            file_choices = "TXT-Datei (*.txt)|*.txt|DAT-Datei (*.dat)|*.dat|Beliebiger Typ (*.*)|*.*"
            filename = self.newname[i] + ".txt"
            dlg = wx.FileDialog(self, message="Daten speichern unter...", defaultFile=filename, wildcard=file_choices, style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
            
            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                self.savename = os.path.splitext(dlg.GetFilename())[0]
                f = open(path, 'w')
                
                if self.cb_header.IsChecked():
                    f.write("----- Header information -------------------------------------------------------\n")
                    f.write("Current name: " + dlg.GetFilename() + "\n")
                    f.write("Original name: " + self.oldname[i] + "\n")
                    f.write("Comment: " + self.comment[i] + "\n")
                    f.write("Original date: " + self.date[i] + "\n")
                    f.write("Wavelength: " + self.wavelength[i] + "\n")
                    f.write("Omega: " + self.omega[i] + "\n")
                    f.write("2Theta: " + self.twotheta[i] + "\n")
                    f.write("Scan type: " + self.scantype[i] + "\n")
                    f.write("Scan axis: " + self.scanaxis[i] + "\n")
                    f.write("First angle: " + str(self.first[i]) + "\n")
                    f.write("Scan range: " + str(self.range[i]) + "\n")
                    f.write("Step width: " + str(self.step[i]) + "\n")
                    f.write("Time per step: " + str(self.time[i]) + "\n")
                    f.write("Number of points: " + str(self.number[i]) + "\n")
                    f.write("--------------------------------------------------------------------------------\n\n")
                   
                f.write("        Angle\t    Intensity\n")
                
                if self.cb_cps.IsChecked():
                    f.write("          deg\t          cps\n")
                    savetxt(f, vstack((self.data[i][:,0], self.data[i][:,1]/self.time[i])).T, fmt='%13.4f', delimiter='\t')
                else:
                    f.write("          deg\t       counts\n")
                    savetxt(f, self.data[i], fmt='%13.4f', delimiter='\t')
                
                f.close()
                self.flash_status_message("Daten gespeichert unter %s" % path)
                
            i = self.list.GetNextSelected(i)
            
    def on_cb_showlog(self, event):
        if self.cb_showlog.IsChecked():
            self.log.Show()
        else:
            self.log.Hide()
        self.vbox_g.Layout()
        
    def on_about(self, event):
        msg = """Anzeige, Konvertierung und Verbinden von Diffraktometer-Dateien:
        
        - Laden einer oder mehrerer njc-, raw-, udf-, x00-, txt- oder dat-Dateien
        - Datenanzeige per Grafik
        - Daten und Grafik speichern
        - Counts pro Sekunde oder Counts
        - Logarithmische y-Achse
        - Kopfdaten mitspeichern
        - Farbe wechseln
        - Daten verbinden
            
        (basiert auf wxPython und matplotlib)
        
        Version 0.5 - 18.07.2013
        """
        dlg = wx.MessageDialog(self, msg, "About", wx.OK)
        dlg.ShowModal()
        dlg.Destroy()
    
    def flash_status_message(self, msg, flash_len_ms=5000):
        self.statusbar.SetStatusText(msg)
        self.timeroff = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.on_flash_status_off, self.timeroff)
        self.timeroff.Start(flash_len_ms, oneShot=True)
        print msg
    
    def on_flash_status_off(self, event):
        self.statusbar.SetStatusText('')

    def on_exit(self, event):
        self.Destroy()
        
        
class RedirectText(object):
    """ Redirects STDOUT and STDERR to log window. """
    def __init__(self, aWxTextCtrl, color):
        self.out = aWxTextCtrl
        self.color = color
 
    def write(self, string):
        self.out.SetDefaultStyle(wx.TextAttr(self.color))
        self.out.WriteText(string)
        self.out.SetInsertionPoint(self.out.GetLastPosition())

        
if __name__ == '__main__':
    app = wx.PySimpleApp()
    app.frame = MainFrame()
    app.frame.Show()
    app.MainLoop()