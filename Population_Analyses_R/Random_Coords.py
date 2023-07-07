#!/usr/bin/env python
# coding: utf-8
# Aidan Shands
import numpy as np
import random
from shapely.geometry import Polygon, Point

# Modified from https://medium.com/the-data-journal/a-quick-trick-to-create-random-lat-long-coordinates-in-python-within-a-defined-polygon-e8997f05123a
# http://apps.headwallphotonics.com/

# Defining areas
SanDiego = Polygon([(33.324980698038814, -117.44899568057035),
(32.85555507188751, -117.1441250751016),
(33.12279080259022, -116.64562043642972),
(33.48662476057229, -116.94774445986722)])

SantaBarbara = Polygon([(34.49596993698643, -119.94461956571561),
(34.41783780990613, -119.88694134305936),
(34.42151958673883, -119.79321423124296),
(34.40027633647306, -119.71047344755155),
(34.42661716394299, -119.60987988065702),
(34.3827111766243, -119.47083416532499),
(34.52765486787859, -119.5429319436453)])

Ventura = Polygon([(34.29738338587903, -119.32820289387854),
(34.46908547354529, -119.33300941243323),
(34.44926989119972, -118.86677711262854),
(34.159427551261174, -118.8976761604801),
(34.12817163164346, -119.14967506184729),
(34.17363096437636, -119.23344581380042),
(34.216794661274804, -119.25129859700354)])

Orange = Polygon([(33.84656390729292, -117.86908177109268),
(33.76839968279126, -117.8697684166005),
(33.771824428378395, -117.7269461509755),
(33.84713418652317, -117.72900608749893)])

Riverside = Polygon([(33.99164363453845, -117.4150212294714),
(33.45595622170995, -117.41639452048703),
(33.467412966961575, -116.71326952048703),
(33.980256653316104, -116.7393620497839)])

ARR = Polygon([(19.23667224763171, -101.74895808076172),
(19.17345014697106, -101.7479281125),
(19.17118018407183, -101.65179774140626),
(19.23667224763171, -101.65488764619141)])

Periban = Polygon([(19.55402665720231, -102.48030735651538),
(19.48736735814167, -102.47687412897632),
(19.489309279692087, -102.3546512285857),
(19.557261845401513, -102.36701084772632)])

SER = Polygon([(19.436863848842464, -101.71316896933594),
(19.359143750660763, -101.71179567832031),
(19.359143750660763, -101.55867373007813),
(19.436863848842464, -101.56897341269531)])

SJNR = Polygon([(19.443979858564077, -102.16516167421118),
(19.443979858564077, -102.10988671083227),
(19.381161736640028, -102.10782677430883),
(19.380837869040377, -102.16413170594946)])

TANR = Polygon([(19.342696490608056, -102.36861869422607),
(19.333058881284224, -102.36801787940674),
(19.332734918147427, -102.35565826026611),
(19.342615504664874, -102.35582992164306)])

TGR = Polygon([(19.760590311773793, -102.52145688474373),
(19.76091341713882, -102.45176236570076),
(19.71405630992815, -102.44935910642342),
(19.713733109743586, -102.52042691648201)])

UPNR_CAV = Polygon([(19.469262218359077, -102.10209861825143),
(19.338438110651563, -102.09042564461862),
(19.347508518642098, -101.97850242684518),
(19.47444120293369, -101.99978843758737)])

ZCA = Polygon([(19.440774845945707, -101.94966331551706),
(19.382489653597126, -101.94760337899362),
(19.385728267763565, -101.8748189551655),
(19.444012299284392, -101.8864919287983)])

ZITR = Polygon([(19.453792536956808, -100.40324804003905),
(19.37737603422256, -100.38058873828123),
(19.41235075600708, -100.28033849414061),
(19.47645153508475, -100.31604406054686)])

# Create function to generate random points within these boundaries
def polygon_random_points (poly, num_points):
    min_x, min_y, max_x, max_y = poly.bounds
    points = []
    while len(points) < num_points:
        random_point = Point([random.uniform(min_x, max_x), random.uniform(min_y, max_y)])
        if (random_point.within(poly)):
            points.append(random_point)
    return points# Choose the number of points desired. This example uses 20 points.

SD_Coords = polygon_random_points(SanDiego,31)
for p in SD_Coords:
    print(p.x,"\t",p.y)

Riverside_Coords = polygon_random_points(Riverside,14)
for p in Riverside_Coords:
    print(p.x,"\t",p.y)

Orange_Coords = polygon_random_points(Orange,3)
for p in Orange_Coords:
    print(p.x,"\t",p.y)

SB_Coords = polygon_random_points(SantaBarbara,10)
for p in SB_Coords:
    print(p.x,"\t",p.y)

Ventura_Coords = polygon_random_points(Ventura,21)
for p in Ventura_Coords:
    print(p.x,"\t",p.y)

ARR_Coords = polygon_random_points(ARR,5)
for p in ARR_Coords:
    print(p.x,"\t",p.y)

Periban_Coords = polygon_random_points(Periban,5)
for p in Periban_Coords:
    print(p.x,"\t",p.y)

SER_Coords = polygon_random_points(SER,5)
for p in SER_Coords:
    print(p.x,"\t",p.y)

SJNR_Coords = polygon_random_points(SJNR,6)
for p in SJNR_Coords:
    print(p.x,"\t",p.y)

TANR_Coords = polygon_random_points(TANR,5)
for p in TANR_Coords:
    print(p.x,"\t",p.y)

TGR_Coords = polygon_random_points(TGR,5)
for p in TGR_Coords:
    print(p.x,"\t",p.y)

UPNR_CAV_Coords = polygon_random_points(UPNR_CAV,9)
for p in UPNR_CAV_Coords:
    print(p.x,"\t",p.y)

ZCA_Coords = polygon_random_points(ZCA,5)
for p in ZCA_Coords:
    print(p.x,"\t",p.y)

ZITR_Coords = polygon_random_points(ZITR,5)
for p in ZITR_Coords:
    print(p.x,"\t",p.y)
