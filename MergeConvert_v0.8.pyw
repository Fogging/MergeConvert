﻿#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2020 by Hartmut Stoecker
# Contact: hartmut.stoecker@physik.tu-freiberg.de
#
# MergeConvert provides a graphical user interface to display X-ray diffraction data and some special tools for X-ray reflectivity measurements.
#
# The present version is 0.8.1, tested on Python 3.8 and wxPython 4.1.

import wx, codecs
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle
 
from pylab import append, arange, argwhere, array, cos, exp, float_, int32, float32, float64, fromfile, linspace, loadtxt, log, mod, ones, pi, savetxt, sin, sum, vstack, zeros
from os.path import splitext, getmtime
from time import gmtime, strftime
from scipy.optimize import leastsq
from copy import deepcopy
from base64 import b64decode
from zlib import decompress, MAX_WBITS
from zipfile import ZipFile


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
        
        nametext = app.frame.filename[self.a] + ' + ' + app.frame.filename[self.b]
        self.namelabel = wx.StaticText(self.panel, -1, "Name:")
        self.name = wx.TextCtrl(self.panel, -1, nametext, style=wx.TE_PROCESS_ENTER, size=(240, 20))  
        self.name.Bind(wx.EVT_TEXT_ENTER, self.do_merge)
        
        self.okButton = wx.Button(self.panel, wx.ID_OK, 'OK', size=(90, 25))
        self.Bind(wx.EVT_BUTTON, self.onOK, self.okButton)
        self.SetAffirmativeId(wx.ID_OK)

        self.closeButton = wx.Button(self.panel, wx.ID_CANCEL, 'Abbrechen', size=(90, 25))
        self.Bind(wx.EVT_BUTTON, self.onClose, self.closeButton)
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
            wx.MessageBox('Verbinden dieser Daten nicht erfolgreich!', 'Fehler', style=wx.ICON_ERROR)
        
        if success:
            self.angle = angle
            self.int = int
            
            if append:
                self.append_data()
            else:
                self.replace_data()
            
    def append_data(self):
        filename = self.name.GetValue()
        date = strftime("%d-%b-%Y, %H:%M:%S")
        comment = app.frame.filename[self.a] + ' + ' + app.frame.filename[self.b] + ' (scaled to ' + str(self.value) + ')'
        wavelength = app.frame.wavelength[self.a]
        radius = app.frame.radius[self.a]
        omega = app.frame.omega[self.a]
        twotheta = app.frame.twotheta[self.a]
        scantype = app.frame.scantype[self.a]
        scanaxis = app.frame.scanaxis[self.a]
        first = self.angle[0]
        range = self.angle[-1] - self.angle[0]
        step = self.angle[1] - self.angle[0]
        time = min(app.frame.time[self.a], app.frame.time[self.b])
        points = len(self.angle)
        data = vstack((self.angle, self.int)).T
        
        app.frame.add_scan(1, filename, date, comment, wavelength, radius, omega, twotheta, scantype, scanaxis, first, range, step, time, points, data)
        
    def replace_data(self):
        i = len(app.frame.name) - 1
        
        app.frame.filename[i] = self.name.GetValue()
        app.frame.name[i] = splitext(self.name.GetValue())[0]
        app.frame.comment[i] = app.frame.filename[self.a] + ' + ' + app.frame.filename[self.b] + ' (scaled to ' + str(self.value) + ')'
        app.frame.date[i] = strftime("%d-%b-%Y, %H:%M:%S")
        app.frame.first[i] = self.angle[0]
        app.frame.range[i] = self.angle[-1] - self.angle[0]
        app.frame.step[i] = self.angle[1] - self.angle[0]
        app.frame.points[i] = len(self.angle)
        app.frame.data[i] = vstack((self.angle, self.int)).T
        
        app.frame.update_list()
        app.frame.draw_figure()
        
    def onOK(self, event):
        self.do_merge(event)
        self.SetReturnCode(wx.ID_OK)
        self.Destroy()
        
    def onClose(self, event):
        app.frame.remove_file(len(app.frame.name)-1)
        app.frame.update_list()
        app.frame.draw_figure()
        
        self.SetReturnCode(wx.ID_CANCEL)
        self.Destroy()


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
        
        nametext = splitext(app.frame.filename[self.a])[0] + '_corr' + splitext(app.frame.filename[self.a])[1]
        self.namelabel = wx.StaticText(self.panel, -1, "Name:")
        self.name = wx.TextCtrl(self.panel, -1, nametext, style=wx.TE_PROCESS_ENTER, size=(260, 20))
        self.name.Bind(wx.EVT_TEXT_ENTER, self.do_correct)
        
        self.okButton = wx.Button(self.panel, wx.ID_OK, 'OK', size=(90, 25))
        self.Bind(wx.EVT_BUTTON, self.onOK, self.okButton)
        self.SetAffirmativeId(wx.ID_OK)

        self.closeButton = wx.Button(self.panel, wx.ID_CANCEL, 'Abbrechen', size=(90, 25))
        self.Bind(wx.EVT_BUTTON, self.onClose, self.closeButton)
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
        
    def append_data(self):
        filename = self.name.GetValue()
        date = strftime("%d-%b-%Y, %H:%M:%S")
        comment = app.frame.filename[self.a] + ' - Corrected using beam size ' + str(self.beamsize) + ' mm and sample size ' + str(self.samplesize) + ' mm'
        wavelength = app.frame.wavelength[self.a]
        radius = app.frame.radius[self.a]
        omega = app.frame.omega[self.a]
        twotheta = app.frame.twotheta[self.a]
        scantype = app.frame.scantype[self.a]
        scanaxis = app.frame.scanaxis[self.a]
        first = self.angle[0]
        range = self.angle[-1] - self.angle[0]
        step = self.angle[1] - self.angle[0]
        time = app.frame.time[self.a]
        points = len(self.angle)
        data = vstack((self.angle, self.res)).T
        
        app.frame.add_scan(1, filename, date, comment, wavelength, radius, omega, twotheta, scantype, scanaxis, first, range, step, time, points, data)
        
    def replace_data(self):
        i = len(app.frame.name) - 1
        
        app.frame.filename[i] = self.name.GetValue()
        app.frame.name[i] = splitext(self.name.GetValue())[0]
        app.frame.comment[i] = app.frame.filename[self.a] + ' - Corrected using beam size ' + str(self.beamsize) + ' mm and sample size ' + str(self.samplesize) + ' mm'
        app.frame.date[i] = strftime("%d-%b-%Y, %H:%M:%S")
        app.frame.data[i] = vstack((self.angle, self.res)).T
        
        app.frame.update_list()
        app.frame.draw_figure()
        
    def onOK(self, event):
        self.do_correct(event)
        self.SetReturnCode(wx.ID_OK)
        self.Destroy()
        
    def onClose(self, event):
        app.frame.remove_file(len(app.frame.name)-1)
        app.frame.update_list()
        app.frame.draw_figure()
        
        self.SetReturnCode(wx.ID_CANCEL)
        self.Destroy()

        
class ScaleWindow(wx.Dialog):
    """ Dialogue to scale and shift datasets. """

    def __init__(self, parent, id, title):
        self.num = app.frame.list.GetSelectedItemCount()
        
        self.selected = []
        i = app.frame.list.GetFirstSelected()
        for j in arange(self.num):
            self.selected.append(i)
            i = app.frame.list.GetNextSelected(i)
            
        self.xshift = deepcopy(app.frame.xshift)
        self.yshift = deepcopy(app.frame.yshift)
        self.scale = deepcopy(app.frame.scale)
        self.displace = deepcopy(app.frame.displace)
        self.radius = deepcopy(app.frame.radius)
            
        height = 140 + 32 * self.num
        
        wx.Dialog.__init__(self, parent, id, title, size=(950, height))
        
        vbox = wx.BoxSizer(wx.VERTICAL)
        box = wx.StaticBox(self, -1, 'Parameter angeben', (5, 5), (930, height-85))
        bsizer = wx.StaticBoxSizer(box, wx.VERTICAL)
        
        t1 = wx.StaticText(self, -1, ' ', size=(370,-1))
        t2 = wx.StaticText(self, -1, ' x-Verschiebung', size=(100,-1))
        t3 = wx.StaticText(self, -1, ' y-Verschiebung', size=(100,-1))
        t4 = wx.StaticText(self, -1, ' Skalierung', size=(100,-1))
        t5 = wx.StaticText(self, -1, u' Höhenfehler (mm)', size=(100,-1))
        t6 = wx.StaticText(self, -1, ' Radius (mm)', size=(100,-1))
        
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        hbox1.Add(t1, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        hbox1.Add(t2, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        hbox1.Add(t3, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        hbox1.Add(t4, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        hbox1.Add(t5, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        hbox1.Add(t6, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        
        bsizer.Add(hbox1, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL)
        
        for i in self.selected:
            name = wx.StaticText(self, -1, app.frame.filename[i] + ':', size=(370,-1))
            t_xshift = wx.TextCtrl(self, 10*i+1, str(app.frame.xshift[i]), style=wx.TE_PROCESS_ENTER, size=(100,20))
            t_yshift = wx.TextCtrl(self, 10*i+2, str(app.frame.yshift[i]), style=wx.TE_PROCESS_ENTER, size=(100,20))
            t_scale = wx.TextCtrl(self, 10*i+3, str(app.frame.scale[i]), style=wx.TE_PROCESS_ENTER, size=(100,20))
            t_displace = wx.TextCtrl(self, 10*i+4, str(app.frame.displace[i]), style=wx.TE_PROCESS_ENTER, size=(100,20))
            t_radius = wx.TextCtrl(self, 10*i+5, str(app.frame.radius[i]), style=wx.TE_PROCESS_ENTER, size=(100,20))
            
            t_xshift.Bind(wx.EVT_TEXT_ENTER, self.refresh)
            t_yshift.Bind(wx.EVT_TEXT_ENTER, self.refresh)
            t_scale.Bind(wx.EVT_TEXT_ENTER, self.refresh)
            t_displace.Bind(wx.EVT_TEXT_ENTER, self.refresh)
            t_radius.Bind(wx.EVT_TEXT_ENTER, self.refresh)

            hbox1 = wx.BoxSizer(wx.HORIZONTAL)
            hbox1.Add(name, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
            hbox1.Add(t_xshift, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
            hbox1.Add(t_yshift, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
            hbox1.Add(t_scale, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
            hbox1.Add(t_displace, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
            hbox1.Add(t_radius, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
            
            bsizer.Add(hbox1, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL)
        
        bsizer.AddSpacer(5)

        okButton = wx.Button(self, wx.ID_OK, 'OK', size=(100, 25))
        self.Bind(wx.EVT_BUTTON, self.onOK, okButton)
        refreshButton = wx.Button(self, -1, 'Aktualisieren', size=(100, 25))
        self.Bind(wx.EVT_BUTTON, self.refresh, refreshButton)
        closeButton = wx.Button(self, wx.ID_CANCEL, 'Abbrechen', size=(100, 25))
        self.Bind(wx.EVT_BUTTON, self.onClose, closeButton)
        
        self.SetAffirmativeId(wx.ID_OK)
        self.SetEscapeId(wx.ID_CANCEL)
        
        hbox3 = wx.BoxSizer(wx.HORIZONTAL)
        hbox3.Add(okButton, 1)
        hbox3.Add(refreshButton, 1, wx.LEFT, 5)
        hbox3.Add(closeButton, 1, wx.LEFT, 5)

        vbox.Add(bsizer, 0, border=7, flag = wx.ALIGN_CENTER | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        vbox.Add(hbox3, 1, border=7, flag = wx.ALIGN_CENTER | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        self.SetSizer(vbox)
        
    def onOK(self, event):
        self.refresh(event)
        self.SetReturnCode(wx.ID_OK)
        self.Destroy()
        
    def refresh(self, event):
        for i in self.selected:
            raw_value = self.FindWindowById(10*i+1).GetValue().strip()
            try:
                app.frame.xshift[i] = float(raw_value)
            except ValueError:
                self.FindWindowById(10*i+1).ChangeValue(str(app.frame.xshift[i]))
        
            raw_value = self.FindWindowById(10*i+2).GetValue().strip()
            try:
                app.frame.yshift[i] = float(raw_value)
            except ValueError:
                self.FindWindowById(10*i+2).ChangeValue(str(app.frame.yshift[i]))
        
            raw_value = self.FindWindowById(10*i+3).GetValue().strip()
            try:
                app.frame.scale[i] = float(raw_value)
            except ValueError:
                self.FindWindowById(10*i+3).ChangeValue(str(app.frame.scale[i]))
        
            raw_value = self.FindWindowById(10*i+4).GetValue().strip()
            try:
                app.frame.displace[i] = float(raw_value)
            except ValueError:
                self.FindWindowById(10*i+4).ChangeValue(str(app.frame.displace[i]))
        
            raw_value = self.FindWindowById(10*i+5).GetValue().strip()
            try:
                app.frame.radius[i] = float(raw_value)
            except ValueError:
                self.FindWindowById(10*i+5).ChangeValue(str(app.frame.radius[i]))
        
        app.frame.update_list()
        app.frame.draw_figure()
        
    def onClose(self, event):
        app.frame.xshift = deepcopy(self.xshift)
        app.frame.yshift = deepcopy(self.yshift)
        app.frame.scale = deepcopy(self.scale)
        app.frame.displace = deepcopy(self.displace)
        app.frame.radius = deepcopy(self.radius)
        
        app.frame.update_list()
        app.frame.draw_figure()
        
        self.SetReturnCode(wx.ID_CANCEL)
        self.Destroy()
            
        
class MainFrame(wx.Frame):
    """ The main frame of the application."""
    
    def __init__(self):
        wx.Frame.__init__(self, None, -1, 'MergeConvert - Anzeige, Konvertierung und Verbinden von Diffraktometer-Dateien')
        
        self.filename = []
        self.name = []
        self.color = []
        self.date = []
        self.comment = []
        self.wavelength = []
        self.omega = []
        self.twotheta = []
        self.scantype = []
        self.scanaxis = []
        self.first = []
        self.range = []
        self.step = []
        self.time = []
        self.points = []
        self.data = []
        self.xshift = []
        self.yshift = []
        self.scale = []
        self.displace = []
        self.radius = []
        self.scaling = []
        
        self.savename = 'Plot'
        self.checked = []
        self.lastcolor = 0
        self.dblclick = 0
        self.press = None
        self.defaultcolors = array([(0,0,1), (0,0.5,0), (1,0,0), (0,0.75,0.75), (0.75,0,0.75), (0.75,0.75,0), (0,0,0), (0.5,0.5,0.5), (0.5,0.75,0.25), (1.0,0.5,0.0), (0.4,0,0.4)])
        
        self.create_menu()
        self.create_status_bar()
        self.create_main_panel()
        self.draw_figure()

    def create_menu(self):
        self.menubar = wx.MenuBar()
        
        menu_file = wx.Menu()
        m_load = menu_file.Append(wx.ID_OPEN, "Datei &laden\tStrg+O", "Daten aus Datei laden")
        self.Bind(wx.EVT_MENU, self.on_open_file, m_load)
        m_delete = menu_file.Append(wx.ID_DELETE, "Datei &entfernen\tEntf", "Datensatz entfernen")
        self.Bind(wx.EVT_MENU, self.on_delete_file, m_delete)
        m_color = menu_file.Append(-1, "&Farbe wechseln\tStrg+F", "Farbe eines Datensatzes wechseln")
        self.Bind(wx.EVT_MENU, self.on_color, m_color)
        menu_file.AppendSeparator()
        m_scale = menu_file.Append(-1, "Daten &skalieren\tStrg+C", "Datensatz skalieren")
        self.Bind(wx.EVT_MENU, self.on_scale, m_scale)
        m_merge = menu_file.Append(-1, "Daten &verbinden\tStrg+V", u"Zwei Datensätze verbinden")
        self.Bind(wx.EVT_MENU, self.on_merge, m_merge)
        m_subtr = menu_file.Append(-1, "Daten &abziehen\tStrg+A", u"Zwei Datensätze subtrahieren")
        self.Bind(wx.EVT_MENU, self.on_subtr, m_subtr)
        m_correct = menu_file.Append(-1, "&XRR-Korrektur\tStrg+K", "XRR-Daten korrigieren")
        self.Bind(wx.EVT_MENU, self.on_correct, m_correct)
        menu_file.AppendSeparator()
        m_savetext = menu_file.Append(wx.ID_SAVE, "&Daten speichern\tStrg+S", "Datensatz als Textdatei speichern")
        self.Bind(wx.EVT_MENU, self.on_save_text, m_savetext)
        m_saveplot = menu_file.Append(-1, "&Grafik speichern\tStrg+G", "Grafik als Datei speichern")
        self.Bind(wx.EVT_MENU, self.on_save_plot, m_saveplot)
        menu_file.AppendSeparator()
        m_exit = menu_file.Append(wx.ID_CANCEL, "&Beenden\tStrg+X", "Programm verlassen")
        self.Bind(wx.EVT_MENU, self.on_exit, m_exit)
        
        menu_help = wx.Menu()
        m_about = menu_help.Append(wx.ID_HELP, "&Hilfe\tF1", "Hilfe zum Programm")
        self.Bind(wx.EVT_MENU, self.on_about, m_about)
        
        self.accel_tbl = wx.AcceleratorTable([(wx.ACCEL_CTRL, ord('O'), m_load.GetId()),
                                              (wx.ACCEL_NORMAL, wx.WXK_DELETE, m_delete.GetId()),
                                              (wx.ACCEL_CTRL, ord('F'), m_color.GetId()),
                                              (wx.ACCEL_CTRL, ord('C'), m_scale.GetId()),
                                              (wx.ACCEL_CTRL, ord('V'), m_merge.GetId()),
                                              (wx.ACCEL_CTRL, ord('A'), m_subtr.GetId()),
                                              (wx.ACCEL_CTRL, ord('K'), m_correct.GetId()),
                                              (wx.ACCEL_CTRL, ord('S'), m_savetext.GetId()),
                                              (wx.ACCEL_CTRL, ord('G'), m_saveplot.GetId()),
                                              (wx.ACCEL_CTRL, ord('X'), m_exit.GetId()),
                                              (wx.ACCEL_NORMAL, wx.WXK_F1, m_about.GetId()),
                                             ])
        self.SetAcceleratorTable(self.accel_tbl)
        
        self.menubar.Append(menu_file, "&Datei")
        self.menubar.Append(menu_help, "&Hilfe")
        self.SetMenuBar(self.menubar)

    def create_main_panel(self):
        """ Creates the main panel with all the controls on it.
        """
        
        self.spWindow = wx.SplitterWindow(self)
        self.panel = wx.Panel(self.spWindow)
        self.graphpanel = wx.Panel(self.spWindow)
        self.spWindow.SetMinimumPaneSize(100)
        self.spWindow.SplitHorizontally(self.panel, self.graphpanel, 150)
        
        # Create the mpl Figure and FigCanvas objects
        self.dpi = 100
        self.fig = Figure((10, 5), dpi=self.dpi)
        self.fig.set_tight_layout(True)
        self.canvas = FigCanvas(self.graphpanel, -1, self.fig)
        self.canvas.mpl_connect('button_press_event', self.on_Press)
        self.canvas.mpl_connect('motion_notify_event', self.on_Motion)
        self.canvas.mpl_connect('button_release_event', self.on_Release)
        self.canvas.mpl_connect('scroll_event', self.on_Scroll)
        
        self.vbox_g = wx.BoxSizer(wx.VERTICAL)
        self.vbox_g.Add(self.canvas, 1, wx.EXPAND)
        self.graphpanel.SetSizer(self.vbox_g)
        
        # Create the top panel
        self.loadbutton = wx.Button(self.panel, wx.ID_OPEN, "Datei\nladen", size=(70,35))
        self.Bind(wx.EVT_BUTTON, self.on_open_file, self.loadbutton)
        self.loadbutton.SetToolTip(wx.ToolTip("Daten aus Datei laden (Strg+O)"))
              
        self.deletebutton = wx.Button(self.panel, wx.ID_DELETE, "Datei\nentfernen", size=(70,35))
        self.Bind(wx.EVT_BUTTON, self.on_delete_file, self.deletebutton)
        self.deletebutton.SetToolTip(wx.ToolTip("Datensatz entfernen (Entf)"))
        
        self.colorbutton = wx.Button(self.panel, -1, "Farbe\nwechseln", size=(70,35))
        self.Bind(wx.EVT_BUTTON, self.on_color, self.colorbutton)
        self.colorbutton.SetToolTip(wx.ToolTip("Farbe eines Datensatzes wechseln (Strg+F)"))
        
        self.scalebutton = wx.Button(self.panel, -1, "Daten\nskalieren", size=(70,35))
        self.Bind(wx.EVT_BUTTON, self.on_scale, self.scalebutton)
        self.scalebutton.SetToolTip(wx.ToolTip("Datensatz skalieren (Strg+C)"))
        
        self.mergebutton = wx.Button(self.panel, -1, "Daten\nverbinden", size=(70,35))
        self.Bind(wx.EVT_BUTTON, self.on_merge, self.mergebutton)
        self.mergebutton.SetToolTip(wx.ToolTip(u"Zwei Datensätze verbinden (Strg+V)"))
        
        self.subtrbutton = wx.Button(self.panel, -1, "Daten\nabziehen", size=(70,35))
        self.Bind(wx.EVT_BUTTON, self.on_subtr, self.subtrbutton)
        self.subtrbutton.SetToolTip(wx.ToolTip(u"Zwei Datensätze subtrahieren (Strg+A)"))
        
        self.correctbutton = wx.Button(self.panel, -1, "XRR-\nKorrektur", size=(70,35))
        self.Bind(wx.EVT_BUTTON, self.on_correct, self.correctbutton)
        self.correctbutton.SetToolTip(wx.ToolTip("XRR-Daten korrigieren (Strg+K)"))

        self.savetextbutton = wx.Button(self.panel, wx.ID_SAVE, "Daten\nspeichern", size=(70,35))
        self.Bind(wx.EVT_BUTTON, self.on_save_text, self.savetextbutton)
        self.savetextbutton.SetToolTip(wx.ToolTip("Datensatz als Textdatei speichern (Strg+S)"))
        
        self.saveplotbutton = wx.Button(self.panel, -1, "Grafik\nspeichern", size=(70,35))
        self.Bind(wx.EVT_BUTTON, self.on_save_plot, self.saveplotbutton)
        self.saveplotbutton.SetToolTip(wx.ToolTip("Grafik als Datei speichern (Strg+G)"))
        
        self.cb_cps = wx.CheckBox(self.panel, -1, "Counts pro Sekunde", style=wx.ALIGN_LEFT)
        self.cb_cps.SetValue(1)
        self.Bind(wx.EVT_CHECKBOX, self.draw_figure, self.cb_cps)
        
        self.cb_log = wx.CheckBox(self.panel, -1, "y logarithmisch", style=wx.ALIGN_LEFT)
        self.cb_log.SetValue(1)
        self.Bind(wx.EVT_CHECKBOX, self.draw_figure, self.cb_log)
        
        self.cb_grid = wx.CheckBox(self.panel, 101, "Gitternetz anzeigen", style=wx.ALIGN_LEFT)
        self.cb_grid.SetValue(1)
        self.Bind(wx.EVT_CHECKBOX, self.draw_figure, self.cb_grid)
        
        self.cb_legend = wx.CheckBox(self.panel, 102, "Legende anzeigen", style=wx.ALIGN_LEFT)
        self.cb_legend.SetValue(1)
        self.Bind(wx.EVT_CHECKBOX, self.draw_figure, self.cb_legend)
        
        self.list = wx.ListCtrl(self.panel, 103, style=wx.LC_REPORT)
        self.list.EnableCheckBoxes(True)
        self.Bind(wx.EVT_LIST_ITEM_CHECKED, self.draw_figure, self.list)
        self.Bind(wx.EVT_LIST_ITEM_UNCHECKED, self.draw_figure, self.list)
        
        self.list.InsertColumn(0, "Dateiname", width=240)
        self.list.InsertColumn(1, "Startwinkel", width=70)
        self.list.InsertColumn(2, "Endwinkel", width=70)
        self.list.InsertColumn(3, "Schrittweite", width=75)
        self.list.InsertColumn(4, "Messzeit", width=65)
        self.list.InsertColumn(5, "Datum", width=130)
        self.list.InsertColumn(6, "Kommentar", width=200)
        self.list.InsertColumn(7, "Skalierung", width=120)
        
        # Layout with box sizers
        flags = wx.ALIGN_LEFT | wx.ALL
        
        self.vbox1 = wx.BoxSizer(wx.VERTICAL)
        self.vbox1.Add(self.cb_cps, 0, border=3, flag=flags)
        self.vbox1.Add(self.cb_log, 0, border=3, flag=flags)
        
        self.vbox2 = wx.BoxSizer(wx.VERTICAL)
        self.vbox2.Add(self.cb_grid, 0, border=3, flag=flags)
        self.vbox2.Add(self.cb_legend, 0, border=3, flag=flags)
 
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
 
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox.Add(self.loadbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.deletebutton, 0, border=3, flag=flags)
        self.hbox.Add(self.colorbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.scalebutton, 0, border=3, flag=flags)
        self.hbox.Add(self.mergebutton, 0, border=3, flag=flags)
        self.hbox.Add(self.subtrbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.correctbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.savetextbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.saveplotbutton, 0, border=3, flag=flags)
        self.hbox.AddSpacer(20)
        self.hbox.Add(self.vbox1, 0, flag = wx.CENTER)
        self.hbox.Add(self.vbox2, 0, flag = wx.CENTER)
        
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT)
        self.vbox.Add(self.list, 1, border=3, flag = wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND)   
        
        self.panel.SetSizer(self.vbox)
        
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.spWindow, 1, wx.EXPAND)
        self.SetSizer(self.sizer)
        self.sizer.Fit(self)
    
    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()
        self.statusbar.SetFieldsCount(3)
        self.statusbar.SetStatusWidths([-1,100,110])

    def draw_figure(self, event=''):
        """ Redraws the figure. """
        try:
            if event.GetId() == 101 or event.GetId() == 102 or event.GetId() == 103:
                event = 'keep'
        except:
            None
        
        if event == 'keep':
            xl, xm = self.axes.get_xlim()
            yl, ym = self.axes.get_ylim()
        
        self.fig.clear()
        
        self.axes = self.fig.add_subplot(111)
        self.axes.locator_params(axis = 'x', nbins = 15)
        self.axes.grid(self.cb_grid.IsChecked())
        
        checked = 0
        for i in arange(len(self.filename)):
            if self.list.IsItemChecked(i):
                checked += 1
                self.checked[i] = 1
                angle = self.data[i][:,0] - 360 / pi * self.displace[i] / self.radius[i] * cos(self.data[i][:,0] * pi / 360) + self.xshift[i]
                intensity = self.data[i][:,1] * self.scale[i] + self.yshift[i]

                if self.cb_cps.IsChecked() and self.time[i] != 0:
                    self.axes.plot(angle, intensity/self.time[i], label=self.filename[i], color=self.color[i])
                elif self.step[i] == 0:
                    self.axes.plot(angle, intensity, '.', label=self.filename[i], color=self.color[i])
                    self.axes.vlines(angle, [0.1], intensity, color=self.color[i])
                else:
                    self.axes.plot(angle, intensity, label=self.filename[i], color=self.color[i])
            else:
                self.checked[i] = 0
                
        if self.cb_log.IsChecked() and self.data != []:
            self.axes.set_yscale('log')
        
        if self.cb_legend.IsChecked() and checked > 0:
            prop = FontProperties(size=8) 
            self.axes.legend(loc=0, prop=prop)
        
        self.axes.tick_params(axis='both', labelsize=8)
        self.axes.set_xlabel('Winkel (deg)', fontsize=10)
        
        if self.cb_cps.IsChecked():
            self.axes.set_ylabel(u'Intensität (cps)', fontsize=10)
        else:
            self.axes.set_ylabel(u'Intensität (counts)', fontsize=10)
        
        if event == 'keep':
            self.axes.set_xlim(xl, xm)
            self.axes.set_ylim(yl, ym)
            
        self.canvas.draw()
        
    def on_Scroll(self, event):
        x, y = event.xdata, event.ydata
        xl, xm = self.axes.get_xlim()
        yl, ym = self.axes.get_ylim()
        factor = 1 - event.step / 20.0
        # print(x, y, xl, xm, yl, ym, event.step, event.button)
        
        if type(x) != 'float' or type(y) == 'float':
            x = (xl + xm) / 2
            y = (yl + ym) / 2
            factor = 0.7
        
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
        self.dblclick = 0
        
        if hasattr(event, 'dblclick'):
            if event.dblclick:
                self.axes.autoscale()
                self.canvas.draw()
                self.dblclick = 1
                
        if event.button == 1 and event.inaxes == self.axes:
            self.press = event.xdata, event.ydata
            self.rectangle = Rectangle((event.xdata, event.ydata), 0, 0, facecolor="grey", alpha=0.25)
            self.axes.add_patch(self.rectangle)
            
        if event.button == 3 and event.inaxes == self.axes:
            self.press = event.xdata, event.ydata
        
    def on_Motion(self, event):
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
        
        if self.press != None: 
            if event.button == 1 and event.inaxes == self.axes:
                x1, y1 = self.press
                x2 = event.xdata
                y2 = event.ydata
                
                self.rectangle.set_height(y2-y1)
                self.rectangle.set_width(x2-x1)
                self.canvas.draw()
            
            if event.button == 3 and event.inaxes == self.axes:
                x1, y1 = self.press
                x2 = event.xdata
                y2 = event.ydata
                xl, xm = self.axes.get_xlim()
                yl, ym = self.axes.get_ylim()
                        
                self.axes.set_xlim(xl+x1-x2, xm+x1-x2)
                if self.axes.get_yscale() == 'log':
                    self.axes.set_ylim(exp(log(yl)+log(y1)-log(y2)), exp(log(ym)+log(y1)-log(y2)))
                else:
                    self.axes.set_ylim(yl+y1-y2, ym+y1-y2)
                self.canvas.draw()
            
    def on_Release(self, event):
        if event.button == 1 and self.press != None and not self.dblclick:
            self.rectangle.remove()
            
            x1, y1 = self.press
            x2 = event.xdata
            y2 = event.ydata
            
            if x2 > x1:
                self.axes.set_xlim(min(x1,x2), max(x1,x2))
                self.axes.set_ylim(min(y1,y2), max(y1,y2))
            elif x1 > x2:
                self.axes.autoscale()
            
            self.press = None
            self.canvas.draw()
            
    def update_list(self):
        self.list.DeleteAllItems()
        self.scaling = []
        
        self.Unbind(wx.EVT_LIST_ITEM_CHECKED, self.list)
        self.Unbind(wx.EVT_LIST_ITEM_UNCHECKED, self.list)
    
        for i in arange(len(self.filename)):
            self.list.InsertItem(i, self.filename[i])
            self.list.SetItem(i, 1, str(self.first[i]))
            self.list.SetItem(i, 2, str(self.range[i]+self.first[i]))
            self.list.SetItem(i, 3, str(self.step[i]))
            self.list.SetItem(i, 4, str(self.time[i]))
            self.list.SetItem(i, 5, self.date[i])
            self.list.SetItem(i, 6, self.comment[i])
            self.list.SetItemTextColour(i, wx.Colour(int(255*self.color[i][0]), int(255*self.color[i][1]), int(255*self.color[i][2])))
            
            scalecomment = ''
            if self.xshift[i] != 0:
                scalecomment += 'x-shift %s deg' %self.xshift[i]
            if self.yshift[i] != 0:
                if scalecomment != '': scalecomment += ', '
                scalecomment += 'y-shift %s counts' %self.yshift[i]
            if self.scale[i] != 1.0: 
                if scalecomment != '': scalecomment += ', '
                scalecomment += 'scaling %s' %self.scale[i]
            if self.displace[i] != 0:
                if scalecomment != '': scalecomment += ', '
                scalecomment += 'height error %s mm (radius %s mm)' %(self.displace[i], self.radius[i])
            self.list.SetItem(i, 7, scalecomment)
            self.scaling.append(scalecomment)
            
            
            if self.checked[i] == 1:
                self.list.CheckItem(i, True)
            else:
                self.list.CheckItem(i, False)
                
        self.Bind(wx.EVT_LIST_ITEM_CHECKED, self.draw_figure, self.list)
        self.Bind(wx.EVT_LIST_ITEM_UNCHECKED, self.draw_figure, self.list)
    
    def on_color(self, event):
        if self.list.GetSelectedItemCount() == 1:
            i = self.list.GetFirstSelected()
            
            data = wx.ColourData()
            data.SetChooseFull(True)
            data.SetCustomColour(0, (10, 80, 161))
            data.SetCustomColour(1, (46, 144, 40))
            data.SetCustomColour(2, (178, 0, 38))
            data.SetCustomColour(3, (78, 188, 206))
            data.SetCustomColour(4, (230, 110, 1))
            data.SetColour(wx.Colour(int(self.color[i][0]*255), int(self.color[i][1]*255), int(self.color[i][2]*255)))
            dlg = wx.ColourDialog(self, data)
           
            if dlg.ShowModal() == wx.ID_OK:
                res = dlg.GetColourData().Colour
                self.color[i] = (res[0]/255.0, res[1]/255.0, res[2]/255.0)
                self.draw_figure('keep')
                self.update_list()
                
        else:
            wx.MessageBox('Bitte genau einen Datensatz in der Liste markieren!', 'Fehler', style=wx.ICON_ERROR)

    def on_scale(self, event):
        if self.list.GetSelectedItemCount() > 0:
            dlg = ScaleWindow(None, -1, 'Daten skalieren')
            dlg.ShowModal()
        else:
            wx.MessageBox('Bitte mindestens einen Datensatz in der Liste markieren!', 'Fehler', style=wx.ICON_ERROR)
            
    def on_merge(self, event):
        x = self.list.GetSelectedItemCount()
        sel = []
        stepOK = 1
        firstOK = 1
        rangeOK = 1
        
        if x >= 2:
            i = self.list.GetFirstSelected()
            sel.append(i)
            stepwidth = self.step[i]
            firstangle = self.first[i]
            scanrange = self.range[i]
            
            while len(sel) != x:
                i = self.list.GetNextSelected(i)
                sel.append(i)
                if (stepwidth - self.step[i]) > 1e-7:
                    stepOK = 0
                if firstangle != self.first[i]:
                    firstOK = 0
                if scanrange != self.range[i]:
                    rangeOK = 0
        
        if x == 2 and stepOK and not firstOK:
            MergeWindow(None, -1, 'Daten verbinden').ShowModal()
        
        elif x >= 2 and stepOK and firstOK and rangeOK:
            filename = self.filename[sel[0]]
            time = self.time[sel[0]]
            data = 1.0 * self.data[sel[0]]
            
            for i in arange(len(sel)-1)+1:
                filename = filename + ' - ' + self.filename[sel[i]]
                time = time + self.time[sel[i]]
                data[:,1] = data[:,1] + self.data[sel[i]][:,1]
            
            self.add_scan(1, filename, self.date[sel[0]], filename, self.wavelength[sel[0]], self.radius[sel[0]], self.omega[sel[0]], self.twotheta[sel[0]], self.scantype[sel[0]], self.scanaxis[sel[0]], self.first[sel[0]], self.range[sel[0]], self.step[sel[0]], time, self.points[sel[0]], data)
            
        elif x < 2:
            wx.MessageBox('Bitte mindestens zwei Dateien in der Liste markieren!', 'Fehler', style=wx.ICON_ERROR)
            
        else:
            wx.MessageBox(u'Dateien müssen gleiche Schrittweite besitzen!\nZusammenfügen von zwei XRR-Dateien erfordert Überlappungsbereich!', 'Fehler', style=wx.ICON_ERROR)
            
    def on_subtr(self, event):
        x = self.list.GetSelectedItemCount()
        
        if x == 2:
            i = self.list.GetFirstSelected()
            j = self.list.GetNextSelected(i)
            stepOK = (self.step[i] - self.step[j] < 1e-7)
            firstOK = (self.first[i] == self.first[j])
            rangeOK = (self.range[i] == self.range[j])
        
        if x == 2 and stepOK and firstOK and rangeOK:
            filename = self.filename[i] + ' + ' + self.filename[j]
            time = self.time[i] + self.time[j]
            data = 1.0 * self.data[i]
            data[:,1] = self.data[i][:,1] - self.data[j][:,1]
            self.add_scan(1, filename, self.date[i], filename, self.wavelength[i], self.radius[i], self.omega[i], self.twotheta[i], self.scantype[i], self.scanaxis[i], self.first[i], self.range[i], self.step[i], time, self.points[i], data)
            
        elif x < 2 or x > 2:
            wx.MessageBox('Bitte genau zwei Dateien in der Liste markieren!', 'Fehler', style=wx.ICON_ERROR)
            
        else:
            wx.MessageBox(u'Dateien müssen gleichen Messbereich und gleiche Schrittweite besitzen!', 'Fehler', style=wx.ICON_ERROR)
            
    def on_correct(self, event):
        x = self.list.GetSelectedItemCount()
        
        if x == 1:
            CorrectWindow(None, -1, 'XRR-Daten korrigieren').ShowModal()
        else:
            wx.MessageBox('Bitte genau einen Datensatz in der Liste markieren!', 'Fehler', style=wx.ICON_ERROR)
       
    def on_delete_file(self, event):
        x = self.list.GetSelectedItemCount()
        dodelete = 0
        
        if x == 0:
            wx.MessageBox('Bitte mindestens einen Datensatz in der Liste markieren!', 'Fehler', style=wx.ICON_ERROR)
        elif x == 1:
            dodelete = 1
        else:
            dlg = wx.MessageDialog(self, 'Alle ' + str(x) + ' markierten Dateien aus der Liste entfernen?', 'Frage', wx.YES_NO | wx.YES_DEFAULT | wx.ICON_QUESTION)
            if dlg.ShowModal() == wx.ID_YES:
                dodelete = 1
                
        if dodelete:
            for i in arange(len(self.checked)-1, -0.5, -1, dtype=int):
                if self.list.IsSelected(i):
                    self.remove_file(i)
            self.update_list()
            self.draw_figure('keep')
            
    def remove_file(self, i):
        del self.filename[i]
        del self.name[i]
        del self.color[i]
        del self.date[i]
        del self.comment[i]
        del self.wavelength[i]
        del self.omega[i]
        del self.twotheta[i]
        del self.scantype[i]
        del self.scanaxis[i]
        del self.first[i]
        del self.range[i]
        del self.step[i]
        del self.time[i]
        del self.points[i]
        del self.data[i]
        del self.xshift[i]
        del self.yshift[i]
        del self.scale[i]
        del self.displace[i]
        del self.radius[i]
        del self.scaling[i]
        del self.checked[i]
        
    def on_open_file(self, event):
        file_choices = "Daten-Typen|*.raw;*.brml;*.udf;*.x00;*.njc;*.val;*.csv;*.dat;*.txt;*.xy|Bruker RAW Version 3 (*.raw)|*.raw|Bruker BRML (*.brml)|*.brml|Philips (*.udf)|*.udf|Philips (*.x00)|*.x00|Seifert (*.njc)|*.njc|Seifert (*.val)|*.val|Text-Datei (*.dat, *.txt, *.xy)|*.dat;*.txt;*.xy|ICSD Powder Pattern Table (*.csv, *.txt)|*.csv;*.txt|Alle Dateien (*.*)|*.*"
        dlg = wx.FileDialog(self, "Datei laden", "", "", file_choices, wx.FD_OPEN|wx.FD_MULTIPLE)
        
        if dlg.ShowModal() == wx.ID_OK:
            filenames = dlg.GetFilenames()
            paths = dlg.GetPaths()
            
            for i in arange(len(paths)):
                filename = filenames[i]
                path = paths[i]
                ext = splitext(filename)[-1].lstrip('.').lower()
                
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

                    wavelength = str(data64[2])
                    first = round(data64[20], 6)
                    range = round(data64[21], 6) - first
                    step = round(data32[44], 6)
                    points = int(round(range / step)) + 1
                    time = round(data32[46], 6)

                    angle = first + step * arange(points)
                    intens = data32[-40-points:-40] * time
                    data = vstack((angle, intens)).T
                    
                    self.add_scan(numberofscans=1, fname=filename, dat=strftime("%d-%b-%Y, %H:%M:%S", gmtime(getmtime(path))), comm=comment, wl=wavelength, fst=first, rng=range, stp=step, t=time, pts=points, d=data)
                    
                elif ext == 'raw':
                    fd = open(path, 'rb')
                    start = fd.read(7)
                    fd.close()
                    
                    if start == 'RAW1.01':
                        fd = open(path, 'rb')
                        fd.seek(0)
                        data_32 = fromfile(fd, float32)
                        fd.seek(0)
                        data_64 = fromfile(fd, float64)
                        fd.seek(4)
                        data_64_shift = fromfile(fd, float64)
                        fd.seek(0)
                        data_i = fromfile(fd, int32)
                        fd.close()

                        fd = codecs.open(path, encoding='Latin-1')
                        fd.seek(16)
                        data_s1 = fd.read(8)
                        fd.seek(26)
                        data_s2 = fd.read(8)
                        fd.seek(326)
                        data_s3 = fd.read(220)
                        fd.close()
                        
                        date = data_s1 + ', ' + data_s2
                        comment = data_s3.strip(u'\x00')
                        wavelength = str(data_64[77])
                        radius = data_32[141]
                        numberofscans = data_i[3]
                        offset = 712
                        
                        omega = []
                        twotheta = []
                        scantype = []
                        scanaxis = []
                        first = []
                        range = []
                        step = []
                        time = []
                        points = []
                        data = []
                    
                        for i in arange(numberofscans):
                            p = data_i[(offset+4)/4]
                            drive = data_i[(offset+196)/4]
                            supplement = data_i[(offset+256)/4]
                            
                            if mod(offset, 8) == 0:
                                omega.append(str(data_64[(offset+8)/8]))
                                twotheta.append(str(data_64[(offset+16)/8]))
                                s = round(data_64[(offset+176)/8], 6)
                            else:
                                omega.append(str(data_64_shift[(offset+8)/8]))
                                twotheta.append(str(data_64_shift[(offset+16)/8]))
                                s = round(data_64_shift[(offset+176)/8], 6)
                                
                            if drive == 0:
                                scantype.append('Locked Coupled')
                                scanaxis.append('2Theta')
                            elif drive == 1:
                                scantype.append('Unlocked Coupled')
                                scanaxis.append('2Theta')
                            elif drive == 2:
                                scantype.append('Detector Scan')
                                scanaxis.append('2Theta')
                            elif drive == 3:
                                scantype.append('Rocking Curve')
                                scanaxis.append('Omega')
                            elif drive == 9999:
                                scantype.append('Tube Scan')
                                scanaxis.append('Omega')
                            else:
                                scantype.append('???')
                                scanaxis.append('???')
                            
                            if scanaxis[-1] == 'Omega':
                                f = round(float(omega[-1]), 6)
                            else:
                                f = round(float(twotheta[-1]), 6)
                            
                            first.append(f)
                            range.append(s * (p - 1))
                            step.append(s)
                            time.append(round(data_32[(offset+192)/4], 6))
                            points.append(p)
                            
                            angle = f + s * arange(p)
                            index = (offset + 304 + supplement) / 4
                            offset = offset + 304 + supplement + p * 4
                            intens = data_32[index:index+p]
                            data.append(vstack((angle, intens)).T)
                        
                        if numberofscans == 1:
                            self.add_scan(numberofscans, filename, date, comment, wavelength, radius, omega[0], twotheta[0], scantype[0], scanaxis[0], first[0], range[0], step[0], time[0], points[0], data[0])
                        else:
                            self.add_scan(numberofscans, filename, date, comment, wavelength, radius, omega, twotheta, scantype, scanaxis, first, range, step, time, points, data)
                            
                    elif start == 'RAW4.00':
                        fd = open(path, 'rb')
                        fd.seek(3)
                        data_32 = fromfile(file=fd, dtype=float32)
                        fd.seek(3)
                        data_64 = fromfile(file=fd, dtype=float64)
                        fd.seek(3)
                        data_i = fromfile(file=fd, dtype=int32)
                        fd.close()
                        
                        fd = codecs.open(path, encoding='Latin-1')
                        fd.seek(12)
                        data_s1 = fd.read(10)
                        fd.seek(24)
                        data_s2 = fd.read(8)
                        fd.close()
                        
                        date = data_s1 + ', ' + data_s2
                        wavelength = str(data_64[61])
                        first = data_64[318]
                        last = data_64[319]
                        points = data_i[374]
                        time = data_32[396]
                        
                        range = last - first
                        step = range / (points-1)
                        
                        angle = linspace(first, last, points)
                        data = vstack((angle, data_32[699:])).T
                    
                        self.add_scan(numberofscans=1, fname=filename, dat=date, wl=wavelength, fst=first, rng=range, stp=step, t=time, pts=points, d=data)
                    
                    else:
                        wx.MessageBox('Beim Laden der Datei:\n\n' + path + '\n\ntrat folgender Fehler auf:\n\nRAW-Version dieser Datei wird nicht verstanden!', 'Fehler beim Laden der Datei', style=wx.ICON_ERROR)
                    
                elif ext == 'brml':
                    f = open(path, 'rb')
                    x = f.read(5)
                    f.close()
                    
                    if x[:4] == b'PK\x03\x04':
                        z = ZipFile(path).read('Experiment0/RawData0.xml').decode('UTF-8')
                    
                    elif x[:5] == '<?xml':
                        root = ET.parse(path).getroot()
                        for containers in root.iter('Containers'):
                            if containers.find('TypeDesc').text == 'BrukerAXS.Common.ExperimentV5.Experiment.Data.DataContainer':
                                y = b64decode(containers.find('SerializedObject').text)
                                z = decompress(y, MAX_WBITS|16).decode('UTF-8')
                                break
                                               
                    angle = []
                    intens = []
                    
                    root = ET.fromstring(z)
                    for datum in root.iter('Datum'):
                        values = datum.text.split(',')
                        angle.append(values[2])
                        intens.append(values[-1])
                    
                    data = vstack((float_(angle), float_(intens))).T
                    
                    try:
                        wavelength = root.find('.//WaveLengthAverage').get('Value')
                        radius = float(root.find('.//SecondaryTracks').find('.//Radius').get('Value'))
                    except:
                        wavelength = '???'
                        radius = 100.0
                    
                    scantype = z.split('ScanName="')[1].split('"')[0]
                    scanaxis = z.split('AxisName="')[1].split('"')[0]
                    first = float(z.split('<Start>')[1].split('</Start>')[0])
                    range = float(z.split('<Stop>')[1].split('</Stop>')[0]) - first
                    step = float(z.split('<Increment>')[1].split('</Increment>')[0])
                    time = float(z.split('MeasuredTimePerStep="')[1].split('"')[0])
                    points = float(z.split('<MeasurementPoints>')[1].split('</MeasurementPoints>')[0])
                    
                    self.add_scan(numberofscans=1, fname=filename, dat=strftime("%d-%b-%Y, %H:%M:%S", gmtime(getmtime(path))), wl=wavelength, rad=radius, stype=scantype, saxis=scanaxis, fst=first, rng=range, stp=step, t=time, pts=points, d=data)
                
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
                        if data_start:
                            for x in line.rstrip('/\n').split(','):
                                if not x.strip() == '':
                                    intens.append(x.strip())
                        if "RawScan" in line:
                            data_start = 1
                    f.close()
                    
                    wavelength = str((l1 + l2*ratio) / (1 + ratio))
                    points = int(round(range / step)) + 1
                    
                    angle = first + step * arange(points)
                    data = vstack((angle, float_(intens))).T
                    
                    self.add_scan(1, filename, date, comment, wavelength, 100.0, '???', '???', scantype, '???', first, range, step, time, points, data)
                    
                elif ext == 'val':
                    seen_intens = 0
                    fn = open(path, 'r')
                    Intensity = array([])
                    while True:
                        zeile = fn.readline()
                        if not zeile: 
                            break
                        if 'FilePar' in zeile:
                            comment1 = fn.readline().strip()
                            comment2 = fn.readline().strip()
                            continue
                        if 'BerPar01' in zeile:
                            omega = fn.readline()
                            twotheta = fn.readline()
                            first = float(twotheta)
                            stop = float(fn.readline())
                            step = float(fn.readline())
                            fn.readline()
                            time = float(fn.readline())
                            fn.readline()
                            points = float(fn.readline())
                            continue
                        if 'Intens' in zeile:
                            seen_intens = 1
                            continue
                        if seen_intens:
                            Intensity = append(Intensity, float(zeile))
                    fn.close()
                    
                    if comment1 != '' and comment2 != '':
                        comment = comment1 + ' - ' + comment2
                    else:
                        comment = comment1 + comment2
                    
                    range = stop - first
                    angle = linspace(first, stop, points)
                    data = vstack((angle, Intensity * time)).T
                    
                    self.add_scan(numberofscans=1, fname=filename, dat=strftime("%d-%b-%Y, %H:%M:%S", gmtime(getmtime(path))), comm=comment, om=omega, tt=twotheta, fst=first, rng=range, stp=step, t=time, pts=points, d=data)
                    
                elif ext == 'x00':
                    f = open(path, 'r')
                    header_length = 1
                    
                    for line in f:
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
                            points = int(line.split()[1])
                        if "ScanData" in line:
                            break
                        else:
                            header_length += 1
                    f.close()
                    
                    angle = first + step * arange(points)
                    intens = loadtxt(path, skiprows=header_length) * time
                    data = vstack((angle, intens)).T
                    
                    self.add_scan(1, filename, date, comment, wavelength, 100.0, omega, twotheta, scantype, scanaxis, first, range, step, time, points, data)
                
                else:
                    try:
                        use_cps = 0
                        ICSD_Pattern = 0
                        header_length = 0
                        time = 0
                        radius = 100.0
                        
                        f = open(path, 'r')
                        
                        for line in f:
                            splitting = line.lstrip(line.split(': ')[0] + ': ').strip()
                            if "Comment:" in line:
                                comment = splitting
                            if "Original date:" in line:
                                date = splitting
                            if "Wavelength:" in line:
                                wavelength = splitting
                            if "Goniometer radius" in line:
                                radius = float(splitting)
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
                            if "H	K	L	2THETA	D-VALUE	MULT	INTENSITY" in line:
                                ICSD_Pattern = 1
                                
                            try:
                                data = loadtxt(path, skiprows=header_length)
                                break
                            except:
                                header_length += 1
                        
                        f.close()
                        
                        if ICSD_Pattern:
                            data = vstack((data[:,3], data[:,6])).T
                            step = 0
                        else:
                            step = data[1,0] - data[0,0]
                            
                        first = data[0,0]
                        range = data[-1,0] - data[0,0]
                        points = len(data[:,0])
                        
                        if use_cps and time != 0:
                            data[:,1] = data[:,1] * time
                        
                        try:
                            self.add_scan(1, filename, date, comment, wavelength, radius, omega, twotheta, scantype, scanaxis, first, range, step, time, points, data)
                        except:
                            self.add_scan(numberofscans=1, fname=filename, dat=strftime("%d-%b-%Y, %H:%M:%S", gmtime(getmtime(path))), fst=first, rng=range, stp=step, t=time, pts=points, d=data)
                            
                    except Exception as error:
                        wx.MessageBox('Beim Laden der Datei:\n\n' + path + '\n\ntrat folgender Fehler auf:\n\n' + str(error), 'Fehler beim Laden der Datei', style=wx.ICON_ERROR)
                        
    def add_scan(self, numberofscans=0, fname='???', dat='???', comm='???', wl='???', rad=100.0, om='???', tt='???', stype='???', saxis='???', fst=0.0, rng=0.0, stp=0.0, t=0.0, pts=0, d=[[0],[1]]):
        if numberofscans > 1:
            dlg = wx.MessageDialog(self, u'Die Datei enthält ' + str(numberofscans) + u' Datensätze.\nAlle Datensätze summieren?', 'Frage', wx.YES_NO | wx.NO_DEFAULT | wx.ICON_QUESTION)
            if dlg.ShowModal() == wx.ID_YES:
                numberofscans = 1
                
                om = om[0]
                tt = tt[0]
                stype = stype[0]
                saxis = saxis[0]
                fst = fst[0]
                rng = rng[0]
                stp = stp[0]
                t = sum(t)
                pts = pts[0]
                
                sumdata = sum(d, axis=0)
                d = d[0]
                d[:,1] = sumdata[:,1]
                
        if numberofscans > 0:
            name = splitext(fname)[0]
            
            if numberofscans == 1:
                self.filename.append(fname)
                self.name.append(name)
                self.checked.append(1)
                self.xshift.append(0.0)
                self.yshift.append(0.0)
                self.scale.append(1.0)
                self.displace.append(0.0)
                
                color = self.defaultcolors[self.lastcolor]
                self.lastcolor = mod( self.lastcolor + 1, len(self.defaultcolors) )
                self.color.append(color)
                
                self.date.append(dat)
                self.comment.append(comm)
                self.wavelength.append(wl)
                self.radius.append(rad)
                self.omega.append(om)
                self.twotheta.append(tt)
                self.scantype.append(stype)
                self.scanaxis.append(saxis)
                self.first.append(fst)
                self.range.append(rng)
                self.step.append(stp)
                self.time.append(t)
                self.points.append(pts)
                self.data.append(d)
            
            else:
                for i in range(numberofscans):
                    self.filename.append('%s - %03i' %(fname,i+1))
                    self.name.append('%s - %03i' %(name,i+1))
                    self.date.append(dat)
                    self.comment.append(comm)
                    self.wavelength.append(wl)
                    self.radius.append(rad)
                
                self.checked.extend(ones(numberofscans, dtype=int))
                self.xshift.extend(zeros(numberofscans, dtype=float))
                self.yshift.extend(zeros(numberofscans, dtype=float))
                self.scale.extend(ones(numberofscans, dtype=float))
                self.displace.extend(zeros(numberofscans, dtype=float))
                
                color = self.defaultcolors[mod( arange(numberofscans) + self.lastcolor, len(self.defaultcolors) )]
                self.lastcolor = mod( numberofscans + self.lastcolor + 1, len(self.defaultcolors) )
                self.color.extend(color)
                
                self.omega.extend(om)
                self.twotheta.extend(tt)
                self.scantype.extend(stype)
                self.scanaxis.extend(saxis)
                self.first.extend(fst)
                self.range.extend(rng)
                self.step.extend(stp)
                self.time.extend(t)
                self.points.extend(pts)
                self.data.extend(d)
                
            self.flash_status_message("%s geladen." %fname)
            self.update_list()
            self.draw_figure()
        
    def on_save_plot(self, event):
        file_choices = "Portable Network Graphics (*.png)|*.png|Encapsulated Postscript (*.eps)|*.eps|Portable Document Format (*.pdf)|*.pdf|Scalable Vector Graphics (*.svg)|*.svg"
        filename = self.savename + '.png'
        dlg = wx.FileDialog(self, message="Grafik speichern unter...", defaultFile=filename, wildcard=file_choices, style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.savename = splitext(dlg.GetFilename())[0]
            self.canvas.print_figure(path, dpi=3*self.dpi)
            self.flash_status_message("Grafik gespeichert unter %s" % path)
        
    def on_save_text(self, event):
        i = self.list.GetFirstSelected()
        
        if i <= -1:
            wx.MessageBox('Bitte mindestens einen Datensatz in der Liste markieren!', 'Fehler', style=wx.ICON_ERROR)
        
        while i > -1:
            file_choices = "TXT-Datei mit Header (*.txt)|*.txt|DAT-Datei ohne Header (*.dat)|*.dat|XY-Datei nur Werte (*.xy)|*.xy"
            filename = self.name[i] + ".txt"
            dlg = wx.FileDialog(self, message="Daten speichern unter...", defaultFile=filename, wildcard=file_choices, style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
            
            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                filename = dlg.GetFilename()
                self.savename = splitext(filename)[0]
                ext = splitext(filename)[-1].lstrip('.').lower()
                f = open(path, 'w')
                
                if ext == 'txt':
                    f.write("----- Header information -------------------------------------------------------\n")
                    f.write("Current name: " + dlg.GetFilename() + "\n")
                    f.write("Original name: " + self.filename[i] + "\n")
                    
                    if self.scaling[i] == '':
                        f.write("Comment: " + self.comment[i] + "\n")
                    else:
                        f.write("Comment: " + self.comment[i] + " - Modifications: " + self.scaling[i] + "\n")
                    
                    f.write("Original date: " + self.date[i] + "\n")
                    f.write("Wavelength: " + self.wavelength[i] + "\n")
                    f.write("Goniometer radius: " + str(self.radius[i]) + "\n")
                    f.write("Omega: " + self.omega[i] + "\n")
                    f.write("2Theta: " + self.twotheta[i] + "\n")
                    f.write("Scan type: " + self.scantype[i] + "\n")
                    f.write("Scan axis: " + self.scanaxis[i] + "\n")
                    f.write("First angle: " + str(self.first[i]) + "\n")
                    f.write("Scan range: " + str(self.range[i]) + "\n")
                    f.write("Step width: " + str(self.step[i]) + "\n")
                    f.write("Time per step: " + str(self.time[i]) + "\n")
                    f.write("Number of points: " + str(self.points[i]) + "\n")
                    f.write("--------------------------------------------------------------------------------\n\n")
                
                if ext == 'txt' or ext == 'dat':
                    f.write("        Angle\t    Intensity\n")
                
                if self.cb_cps.IsChecked():
                    if ext == 'txt' or ext == 'dat':
                        f.write("          deg\t          cps\n")
                    savetxt(f, vstack((self.data[i][:,0], self.data[i][:,1]/self.time[i])).T, fmt='%13.4f', delimiter='\t')
                else:
                    if ext == 'txt' or ext == 'dat':
                        f.write("          deg\t       counts\n")
                    savetxt(f, self.data[i], fmt='%13.4f', delimiter='\t')
                
                f.close()
                self.flash_status_message("Daten gespeichert unter %s" % path)
                
            i = self.list.GetNextSelected(i)
            
    def on_about(self, event):
        msg = """Anzeige, Konvertierung und Verbinden von Diffraktometer-Dateien:
        
        - Laden einer oder mehrerer Dateien
        - Formate: raw, udf, x00, njc, val, dat, txt
        - Daten und Grafik speichern
        - Daten skalieren
        - Daten verbinden
        - XRR-Korrektur bei kleiner Probe
            
        (basiert auf wxPython und matplotlib)
        
        Version 0.8.1 - 05.11.2020
        """
        dlg = wx.MessageDialog(self, msg, "About", wx.OK)
        dlg.ShowModal()
        dlg.Destroy()
    
    def flash_status_message(self, msg, flash_len_ms=5000):
        self.statusbar.SetStatusText(msg)
        self.timeroff = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.on_flash_status_off, self.timeroff)
        self.timeroff.Start(flash_len_ms, oneShot=True)
    
    def on_flash_status_off(self, event):
        self.statusbar.SetStatusText('')

    def on_exit(self, event):
        self.Destroy()
        
        
if __name__ == '__main__':
    app = wx.App(redirect=True)
    app.frame = MainFrame()
    app.frame.Show()
    app.MainLoop()