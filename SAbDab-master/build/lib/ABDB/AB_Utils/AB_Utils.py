'''
Created on 18 Feb 2013

@author: dunbar

Functions, dictionaries and lists for defining regions of the antibody.
General utility functions for ABDB.


Definitions of regions of the antibody in Chothia numbering:
    o Chothia defintions of the CDRs.
    o Kabat defintions of the CDRs.
    o Contact defintions of the CDRs.

ALL numbering is in the Chothia numbering scheme

@change: added in kabat and contact defintions for CDR regions. Chothia region_annotations are default. chothia_region_annotations = region_annotations

'''
import os

aa1="ACDEFGHIKLMNPQRSTVWY"

positions = ['H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H6A', 'H6B', 'H6C', 'H6D', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13', 'H14', 'H15', 'H16', 'H17', 'H18', 'H19', 'H20', 'H21', 'H22', 'H23', 'H24', 'H25', 'H26', 'H27', 'H28', 'H29', 'H30', 'H31', 'H31A', 'H31B', 'H31C', 'H31D', 'H31E', 'H31F', 'H31G', 'H31H', 'H32', 'H33', 'H34', 'H35', 'H36', 'H37', 'H38', 'H39', 'H40', 'H41', 'H42', 'H43', 'H44', 'H45', 'H46', 'H47', 'H48', 'H49', 'H50', 'H51', 'H52', 'H52A', 'H52B', 'H52C', 'H53', 'H54', 'H55', 'H56', 'H57', 'H58', 'H59', 'H60', 'H61', 'H62', 'H63', 'H64', 'H65', 'H66', 'H67', 'H68', 'H69', 'H70', 'H71', 'H72', 'H73', 'H74', 'H75', 'H76', 'H77', 'H78', 'H79', 'H80', 'H81', 'H82', 'H82A', 'H82B', 'H82C', 'H83', 'H84', 'H85', 'H86', 'H87', 'H88', 'H89', 'H90', 'H91', 'H92', 'H93', 'H94', 'H95', 'H96', 'H97', 'H98', 'H99', 'H100', 'H100A', 'H100B', 'H100C', 'H100D', 'H100E', 'H100F', 'H100G', 'H100H', 'H100I', 'H100J', 'H100K', 'H100L', 'H100M', 'H100N', 'H100O', 'H100P', 'H100Q', 'H100R', 'H100S', 'H100T', 'H101', 'H102', 'H103', 'H104', 'H105', 'H106', 'H107', 'H108', 'H109', 'H110', 'H111', 'H112', 'H113', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8', 'L9', 'L10', 'L11', 'L12', 'L13', 'L14', 'L15', 'L16', 'L17', 'L18', 'L19', 'L20', 'L21', 'L22', 'L23', 'L24', 'L25', 'L26', 'L27', 'L28', 'L29', 'L30', 'L30A', 'L30B', 'L30C', 'L30D', 'L30E', 'L30F', 'L31', 'L32', 'L33', 'L34', 'L35', 'L36', 'L37', 'L38', 'L39', 'L40', 'L41', 'L42', 'L43', 'L44', 'L45', 'L46', 'L47', 'L48', 'L49', 'L50', 'L51', 'L52', 'L53', 'L54', 'L54A', 'L54B', 'L54C', 'L54D', 'L54E', 'L55', 'L56', 'L57', 'L58', 'L59', 'L60', 'L61', 'L62', 'L63', 'L64', 'L65', 'L66', 'L66A', 'L66B', 'L66C', 'L66D', 'L66E', 'L66F', 'L66G', 'L66H', 'L67', 'L68', 'L69', 'L70', 'L71', 'L72', 'L73', 'L74', 'L75', 'L76', 'L77', 'L78', 'L79', 'L80', 'L81', 'L82', 'L83', 'L84', 'L85', 'L86', 'L87', 'L88', 'L89', 'L90', 'L91', 'L92', 'L93', 'L94', 'L95', 'L95A', 'L95B', 'L95C', 'L95D', 'L95E', 'L95F', 'L96', 'L97', 'L98', 'L99', 'L100', 'L101', 'L102', 'L103', 'L104', 'L105', 'L106', 'L106A', 'L107', 'L108', 'L109', 'L110', 'L111']#: Known positions in an antibody in abnum format and chothia numbering
 
positions_tuples = {}#: Known positions in an antibody in L{Bio.PDB} tuple format.
positions_tuples["H"] = [(1, ' '), (2, ' '), (3, ' '), (4, ' '), (5, ' '), (6, ' '), (6, 'A'), (6, 'B'), (6, 'C'), (6, 'D'), (7, ' '), (8, ' '), (9, ' '), (10, ' '), (11, ' '), (12, ' '), (13, ' '), (14, ' '), (15, ' '), (16, ' '), (17, ' '), (18, ' '), (19, ' '), (20, ' '), (21, ' '), (22, ' '), (23, ' '), (24, ' '), (25, ' '), (26, ' '), (27, ' '), (28, ' '), (29, ' '), (30, ' '), (31, ' '), (31, 'A'), (31, 'B'), (31, 'C'), (31, 'D'), (31, 'E'), (31, 'F'), (31, 'G'), (31, 'H'), (32, ' '), (33, ' '), (34, ' '), (35, ' '), (36, ' '), (37, ' '), (38, ' '), (39, ' '), (40, ' '), (41, ' '), (42, ' '), (43, ' '), (44, ' '), (45, ' '), (46, ' '), (47, ' '), (48, ' '), (49, ' '), (50, ' '), (51, ' '), (52, ' '), (52, 'A'), (52, 'B'), (52, 'C'), (53, ' '), (54, ' '), (55, ' '), (56, ' '), (57, ' '), (58, ' '), (59, ' '), (60, ' '), (61, ' '), (62, ' '), (63, ' '), (64, ' '), (65, ' '), (66, ' '), (67, ' '), (68, ' '), (69, ' '), (70, ' '), (71, ' '), (72, ' '), (73, ' '), (74, ' '), (75, ' '), (76, ' '), (77, ' '), (78, ' '), (79, ' '), (80, ' '), (81, ' '), (82, ' '), (82, 'A'), (82, 'B'), (82, 'C'), (83, ' '), (84, ' '), (85, ' '), (86, ' '), (87, ' '), (88, ' '), (89, ' '), (90, ' '), (91, ' '), (92, ' '), (93, ' '), (94, ' '), (95, ' '), (96, ' '), (97, ' '), (98, ' '), (99, ' '), (100, ' '), (100, 'A'), (100, 'B'), (100, 'C'), (100, 'D'), (100, 'E'), (100, 'F'), (100, 'G'), (100, 'H'), (100, 'I'), (100, 'J'), (100, 'K'), (100, 'L'), (100, 'M'), (100, 'N'), (100, 'O'), (100, 'P'), (100, 'Q'), (100, 'R'), (100, 'S'), (100, 'T'), (101, ' '), (102, ' '), (103, ' '), (104, ' '), (105, ' '), (106, ' '), (107, ' '), (108, ' '), (109, ' '), (110, ' '), (111, ' '), (112, ' '), (113, ' ')]
positions_tuples["L"] = [(1, ' '), (2, ' '), (3, ' '), (4, ' '), (5, ' '), (6, ' '), (7, ' '), (8, ' '), (9, ' '), (10, ' '), (11, ' '), (12, ' '), (13, ' '), (14, ' '), (15, ' '), (16, ' '), (17, ' '), (18, ' '), (19, ' '), (20, ' '), (21, ' '), (22, ' '), (23, ' '), (24, ' '), (25, ' '), (26, ' '), (27, ' '), (28, ' '), (29, ' '), (30, ' '), (30, 'A'), (30, 'B'), (30, 'C'), (30, 'D'), (30, 'E'), (30, 'F'), (31, ' '), (32, ' '), (33, ' '), (34, ' '), (35, ' '), (36, ' '), (37, ' '), (38, ' '), (39, ' '), (40, ' '), (41, ' '), (42, ' '), (43, ' '), (44, ' '), (45, ' '), (46, ' '), (47, ' '), (48, ' '), (49, ' '), (50, ' '), (51, ' '), (52, ' '), (53, ' '), (54, ' '), (54, 'A'), (54, 'B'), (54, 'C'), (54, 'D'), (54, 'E'), (55, ' '), (56, ' '), (57, ' '), (58, ' '), (59, ' '), (60, ' '), (61, ' '), (62, ' '), (63, ' '), (64, ' '), (65, ' '), (66, ' '), (66, 'A'), (66, 'B'), (66, 'C'), (66, 'D'), (66, 'E'), (66, 'F'), (66, 'G'), (66, 'H'), (67, ' '), (68, ' '), (69, ' '), (70, ' '), (71, ' '), (72, ' '), (73, ' '), (74, ' '), (75, ' '), (76, ' '), (77, ' '), (78, ' '), (79, ' '), (80, ' '), (81, ' '), (82, ' '), (83, ' '), (84, ' '), (85, ' '), (86, ' '), (87, ' '), (88, ' '), (89, ' '), (90, ' '), (91, ' '), (92, ' '), (93, ' '), (94, ' '), (95, ' '), (95, 'A'), (95, 'B'), (95, 'C'), (95, 'D'), (95, 'E'), (95, 'F'), (96, ' '), (97, ' '), (98, ' '), (99, ' '), (100, ' '), (101, ' '), (102, ' '), (103, ' '), (104, ' '), (105, ' '), (106, ' '), (106, 'A'), (107, ' '), (108, ' '), (109, ' '), (110, ' '), (111, ' ')]

# chothia definition for CDRs
region_annotations={"H1":"fwH1",
"H2":"fwH1",
"H3":"fwH1",
"H4":"fwH1",
"H5":"fwH1",
"H6":"fwH1",
"H6A":"fwH1",
"H6B":"fwH1",
"H6C":"fwH1",
"H6D":"fwH1",
"H7":"fwH1",
"H8":"fwH1",
"H9":"fwH1",
"H10":"fwH1",
"H11":"fwH1",
"H12":"fwH1",
"H13":"fwH1",
"H14":"fwH1",
"H15":"fwH1",
"H16":"fwH1",
"H17":"fwH1",
"H18":"fwH1",
"H19":"fwH1",
"H20":"fwH1",
"H21":"fwH1",
"H22":"fwH1",
"H23":"fwH1",
"H24":"fwH1",
"H25":"fwH1",
"H26":"CDRH1",
"H27":"CDRH1",
"H28":"CDRH1",
"H29":"CDRH1",
"H30":"CDRH1",
"H31":"CDRH1",
"H31A":"CDRH1",
"H31B":"CDRH1",
"H31C":"CDRH1",
"H31D":"CDRH1",
"H31E":"CDRH1",
"H31F":"CDRH1",
"H31G":"CDRH1",
"H31H":"CDRH1",
"H32":"CDRH1",
"H33":"fwH2",
"H34":"fwH2",
"H35":"fwH2",
"H36":"fwH2",
"H37":"fwH2",
"H38":"fwH2",
"H39":"fwH2",
"H40":"fwH2",
"H41":"fwH2",
"H42":"fwH2",
"H43":"fwH2",
"H44":"fwH2",
"H45":"fwH2",
"H46":"fwH2",
"H47":"fwH2",
"H48":"fwH2",
"H49":"fwH2",
"H50":"fwH2",
"H51":"fwH2",
"H52":"CDRH2",
"H52A":"CDRH2",
"H52B":"CDRH2",
"H52C":"CDRH2",
"H53":"CDRH2",
"H54":"CDRH2",
"H55":"CDRH2",
"H56":"CDRH2",
"H57":"afwH3",
"H58":"afwH3",
"H59":"afwH3",
"H60":"afwH3",
"H61":"afwH3",
"H62":"afwH3",
"H63":"afwH3",
"H64":"afwH3",
"H65":"afwH3",
"H66":"afwH3",
"H67":"afwH3",
"H68":"afwH3",
"H69":"afwH3",
"H70":"afwH3",
"H71":"afwH3",
"H72":"afwH3",
"H73":"afwH3",
"H74":"afwH3",
"H75":"afwH3",
"H76":"bfwH3",
"H77":"bfwH3",
"H78":"bfwH3",
"H79":"bfwH3",
"H80":"bfwH3",
"H81":"bfwH3",
"H82":"bfwH3",
"H82A":"bfwH3",
"H82B":"bfwH3",
"H82C":"bfwH3",
"H83":"bfwH3",
"H84":"bfwH3",
"H85":"bfwH3",
"H86":"bfwH3",
"H87":"bfwH3",
"H88":"bfwH3",
"H89":"bfwH3",
"H90":"bfwH3",
"H91":"bfwH3",
"H92":"bfwH3",
"H93":"bfwH3",
"H94":"bfwH3",
"H95":"CDRH3",
"H96":"CDRH3",
"H97":"CDRH3",
"H98":"CDRH3",
"H99":"CDRH3",
"H100":"CDRH3",
"H100A":"CDRH3",
"H100B":"CDRH3",
"H100C":"CDRH3",
"H100D":"CDRH3",
"H100E":"CDRH3",
"H100F":"CDRH3",
"H100G":"CDRH3",
"H100H":"CDRH3",
"H100I":"CDRH3",
"H100J":"CDRH3",
"H100K":"CDRH3",
"H100L":"CDRH3",
"H100M":"CDRH3",
"H100N":"CDRH3",
"H100O":"CDRH3",
"H100P":"CDRH3",
"H100Q":"CDRH3",
"H100R":"CDRH3",
"H100S":"CDRH3",
"H100T":"CDRH3",
"H101":"CDRH3",
"H102":"CDRH3",
"H103":"fwH4",
"H104":"fwH4",
"H105":"fwH4",
"H106":"fwH4",
"H107":"fwH4",
"H108":"fwH4",
"H109":"fwH4",
"H110":"fwH4",
"H111":"fwH4",
"H112":"fwH4",
"H113":"fwH4",
"L1":"fwL1",
"L2":"fwL1",
"L3":"fwL1",
"L4":"fwL1",
"L5":"fwL1",
"L6":"fwL1",
"L7":"fwL1",
"L8":"fwL1",
"L9":"fwL1",
"L10":"fwL1",
"L11":"fwL1",
"L12":"fwL1",
"L13":"fwL1",
"L14":"fwL1",
"L15":"fwL1",
"L16":"fwL1",
"L17":"fwL1",
"L18":"fwL1",
"L19":"fwL1",
"L20":"fwL1",
"L21":"fwL1",
"L22":"fwL1",
"L23":"fwL1",
"L24":"CDRL1",
"L25":"CDRL1",
"L26":"CDRL1",
"L27":"CDRL1",
"L28":"CDRL1",
"L29":"CDRL1",
"L30":"CDRL1",
"L30A":"CDRL1",
"L30B":"CDRL1",
"L30C":"CDRL1",
"L30D":"CDRL1",
"L30E":"CDRL1",
"L30F":"CDRL1",
"L31":"CDRL1",
"L32":"CDRL1",
"L33":"CDRL1",
"L34":"CDRL1",
"L35":"fwL2",
"L36":"fwL2",
"L37":"fwL2",
"L38":"fwL2",
"L39":"fwL2",
"L40":"fwL2",
"L41":"fwL2",
"L42":"fwL2",
"L43":"fwL2",
"L44":"fwL2",
"L45":"fwL2",
"L46":"fwL2",
"L47":"fwL2",
"L48":"fwL2",
"L49":"fwL2",
"L50":"CDRL2",
"L51":"CDRL2",
"L52":"CDRL2",
"L53":"CDRL2",
"L54":"CDRL2",
"L54A":"CDRL2",
"L54B":"CDRL2",
"L54C":"CDRL2",
"L54D":"CDRL2",
"L54E":"CDRL2",
"L55":"CDRL2",
"L56":"CDRL2",
"L57":"afwL3",
"L58":"afwL3",
"L59":"afwL3",
"L60":"afwL3",
"L61":"afwL3",
"L62":"afwL3",
"L63":"afwL3",
"L64":"afwL3",
"L65":"afwL3",
"L66":"afwL3",
"L66A":"afwL3",
"L66B":"afwL3",
"L66C":"afwL3",
"L66D":"afwL3",
"L66E":"afwL3",
"L66F":"afwL3",
"L66G":"afwL3",
"L66H":"afwL3",
"L67":"afwL3",
"L68":"afwL3",
"L69":"afwL3",
"L70":"afwL3",
"L71":"afwL3",
"L72":"afwL3",
"L73":"bfwL3",
"L74":"bfwL3",
"L75":"bfwL3",
"L76":"bfwL3",
"L77":"bfwL3",
"L78":"bfwL3",
"L79":"bfwL3",
"L80":"bfwL3",
"L81":"bfwL3",
"L82":"bfwL3",
"L83":"bfwL3",
"L84":"bfwL3",
"L85":"bfwL3",
"L86":"bfwL3",
"L87":"bfwL3",
"L88":"bfwL3",
"L89":"CDRL3",
"L90":"CDRL3",
"L91":"CDRL3",
"L92":"CDRL3",
"L93":"CDRL3",
"L94":"CDRL3",
"L95":"CDRL3",
"L95A":"CDRL3",
"L95B":"CDRL3",
"L95C":"CDRL3",
"L95D":"CDRL3",
"L95E":"CDRL3",
"L95F":"CDRL3",
"L96":"CDRL3",
"L97":"CDRL3",
"L98":"fwL4",
"L99":"fwL4",
"L100":"fwL4",
"L101":"fwL4",
"L102":"fwL4",
"L103":"fwL4",
"L104":"fwL4",
"L105":"fwL4",
"L106":"fwL4",
"L106A":"fwL4",
"L107":"fwL4",
"L108":"fwL4",
"L109":"fwL4",
"L110":"fwL4",
"L111":"fwL4"}#: abnum format to chothia regions ( framework region 3 divided into two sections a and b)

# Chothia definition for CDRs 
regions = {}#: Positions of the antibody grouped into region of the structure using abnum format annotations. Chothia definitions for the loops.
for p in positions:
    try:
        regions[region_annotations[p]].append(p)
    except KeyError:
        regions[region_annotations[p]] = [p]
         
regions["loops"] = regions["CDRH1"]+regions["CDRH2"] +regions["CDRH3"]+regions["CDRL1"]+regions["CDRL2"]+regions["CDRL3"]
regions["framework"] = regions["fwH1"]+regions["fwH2"]+regions["afwH3"]+regions["bfwH3"]+regions["fwH4"]+regions["fwL1"]+regions["fwL2"]+regions["afwL3"]+regions["bfwL3"]+regions["fwL4"] 
# interface positions from 70% non-redundant set, 4.5A any atom contact cutoff and 5% filter of contact frequency.
regions["interface"] = ["H100", "H100A", "H100B", "H100C", "H100D", "H100E", "H100F", "H100G", "H101", "H102", "H103", "H104", "H105", "H106", "H35", "H37", "H39", "H43", "H44", "H45", "H46", "H47", "H50", "H58", "H59", "H60", "H61", "H62", "H91", "H93", "H95", "H96", "H98", "H99", "L1", "L100", "L32", "L34", "L36", "L38", "L41", "L42", "L43", "L44", "L45", "L46", "L49", "L50", "L55", "L85", "L87", "L89", "L91", "L94", "L95", "L95A", "L95B", "L95C", "L96", "L98", "L99"]

chothia_regions=regions
chothia_region_annotations=region_annotations

# Kabat definition for CDRs (in chothia numbering)
kabat_region_annotations={"H1":"fwH1",
"H2":"fwH1",
"H3":"fwH1",
"H4":"fwH1",
"H5":"fwH1",
"H6":"fwH1",
"H6A":"fwH1",
"H6B":"fwH1",
"H6C":"fwH1",
"H6D":"fwH1",
"H7":"fwH1",
"H8":"fwH1",
"H9":"fwH1",
"H10":"fwH1",
"H11":"fwH1",
"H12":"fwH1",
"H13":"fwH1",
"H14":"fwH1",
"H15":"fwH1",
"H16":"fwH1",
"H17":"fwH1",
"H18":"fwH1",
"H19":"fwH1",
"H20":"fwH1",
"H21":"fwH1",
"H22":"fwH1",
"H23":"fwH1",
"H24":"fwH1",
"H25":"fwH1",
"H26":"fwH1",
"H27":"fwH1",
"H28":"fwH1",
"H29":"fwH1",
"H30":"fwH1",
"H31":"CDRH1",
"H31A":"CDRH1",
"H31B":"CDRH1",
"H31C":"CDRH1",
"H31D":"CDRH1",
"H31E":"CDRH1",
"H31F":"CDRH1",
"H31G":"CDRH1",
"H31H":"CDRH1",
"H32":"CDRH1",
"H33":"CDRH1",
"H34":"CDRH1",
"H35":"CDRH1",
"H36":"fwH2",
"H37":"fwH2",
"H38":"fwH2",
"H39":"fwH2",
"H40":"fwH2",
"H41":"fwH2",
"H42":"fwH2",
"H43":"fwH2",
"H44":"fwH2",
"H45":"fwH2",
"H46":"fwH2",
"H47":"fwH2",
"H48":"fwH2",
"H49":"fwH2",
"H50":"CDRH2",
"H51":"CDRH2",
"H52":"CDRH2",
"H52A":"CDRH2",
"H52B":"CDRH2",
"H52C":"CDRH2",
"H53":"CDRH2",
"H54":"CDRH2",
"H55":"CDRH2",
"H56":"CDRH2",
"H57":"CDRH2",
"H58":"CDRH2",
"H59":"CDRH2",
"H60":"CDRH2",
"H61":"CDRH2",
"H62":"CDRH2",
"H63":"CDRH2",
"H64":"CDRH2",
"H65":"CDRH2",
"H66":"afwH3",
"H67":"afwH3",
"H68":"afwH3",
"H69":"afwH3",
"H70":"afwH3",
"H71":"afwH3",
"H72":"afwH3",
"H73":"afwH3",
"H74":"afwH3",
"H75":"afwH3",
"H76":"bfwH3",
"H77":"bfwH3",
"H78":"bfwH3",
"H79":"bfwH3",
"H80":"bfwH3",
"H81":"bfwH3",
"H82":"bfwH3",
"H82A":"bfwH3",
"H82B":"bfwH3",
"H82C":"bfwH3",
"H83":"bfwH3",
"H84":"bfwH3",
"H85":"bfwH3",
"H86":"bfwH3",
"H87":"bfwH3",
"H88":"bfwH3",
"H89":"bfwH3",
"H90":"bfwH3",
"H91":"bfwH3",
"H92":"bfwH3",
"H93":"bfwH3",
"H94":"bfwH3",
"H95":"CDRH3",
"H96":"CDRH3",
"H97":"CDRH3",
"H98":"CDRH3",
"H99":"CDRH3",
"H100":"CDRH3",
"H100A":"CDRH3",
"H100B":"CDRH3",
"H100C":"CDRH3",
"H100D":"CDRH3",
"H100E":"CDRH3",
"H100F":"CDRH3",
"H100G":"CDRH3",
"H100H":"CDRH3",
"H100I":"CDRH3",
"H100J":"CDRH3",
"H100K":"CDRH3",
"H100L":"CDRH3",
"H100M":"CDRH3",
"H100N":"CDRH3",
"H100O":"CDRH3",
"H100P":"CDRH3",
"H100Q":"CDRH3",
"H100R":"CDRH3",
"H100S":"CDRH3",
"H100T":"CDRH3",
"H101":"CDRH3",
"H102":"CDRH3",
"H103":"fwH4",
"H104":"fwH4",
"H105":"fwH4",
"H106":"fwH4",
"H107":"fwH4",
"H108":"fwH4",
"H109":"fwH4",
"H110":"fwH4",
"H111":"fwH4",
"H112":"fwH4",
"H113":"fwH4",
"L1":"fwL1",
"L2":"fwL1",
"L3":"fwL1",
"L4":"fwL1",
"L5":"fwL1",
"L6":"fwL1",
"L7":"fwL1",
"L8":"fwL1",
"L9":"fwL1",
"L10":"fwL1",
"L11":"fwL1",
"L12":"fwL1",
"L13":"fwL1",
"L14":"fwL1",
"L15":"fwL1",
"L16":"fwL1",
"L17":"fwL1",
"L18":"fwL1",
"L19":"fwL1",
"L20":"fwL1",
"L21":"fwL1",
"L22":"fwL1",
"L23":"fwL1",
"L24":"CDRL1",
"L25":"CDRL1",
"L26":"CDRL1",
"L27":"CDRL1",
"L28":"CDRL1",
"L29":"CDRL1",
"L30":"CDRL1",
"L30A":"CDRL1",
"L30B":"CDRL1",
"L30C":"CDRL1",
"L30D":"CDRL1",
"L30E":"CDRL1",
"L30F":"CDRL1",
"L31":"CDRL1",
"L32":"CDRL1",
"L33":"CDRL1",
"L34":"CDRL1",
"L35":"fwL2",
"L36":"fwL2",
"L37":"fwL2",
"L38":"fwL2",
"L39":"fwL2",
"L40":"fwL2",
"L41":"fwL2",
"L42":"fwL2",
"L43":"fwL2",
"L44":"fwL2",
"L45":"fwL2",
"L46":"fwL2",
"L47":"fwL2",
"L48":"fwL2",
"L49":"fwL2",
"L50":"CDRL2",
"L51":"CDRL2",
"L52":"CDRL2",
"L53":"CDRL2",
"L54":"CDRL2",
"L54A":"CDRL2",
"L54B":"CDRL2",
"L54C":"CDRL2",
"L54D":"CDRL2",
"L54E":"CDRL2",
"L55":"CDRL2",
"L56":"CDRL2",
"L57":"afwL3",
"L58":"afwL3",
"L59":"afwL3",
"L60":"afwL3",
"L61":"afwL3",
"L62":"afwL3",
"L63":"afwL3",
"L64":"afwL3",
"L65":"afwL3",
"L66":"afwL3",
"L66A":"afwL3",
"L66B":"afwL3",
"L66C":"afwL3",
"L66D":"afwL3",
"L66E":"afwL3",
"L66F":"afwL3",
"L66G":"afwL3",
"L66H":"afwL3",
"L67":"afwL3",
"L68":"afwL3",
"L69":"afwL3",
"L70":"afwL3",
"L71":"afwL3",
"L72":"afwL3",
"L73":"bfwL3",
"L74":"bfwL3",
"L75":"bfwL3",
"L76":"bfwL3",
"L77":"bfwL3",
"L78":"bfwL3",
"L79":"bfwL3",
"L80":"bfwL3",
"L81":"bfwL3",
"L82":"bfwL3",
"L83":"bfwL3",
"L84":"bfwL3",
"L85":"bfwL3",
"L86":"bfwL3",
"L87":"bfwL3",
"L88":"bfwL3",
"L89":"CDRL3",
"L90":"CDRL3",
"L91":"CDRL3",
"L92":"CDRL3",
"L93":"CDRL3",
"L94":"CDRL3",
"L95":"CDRL3",
"L95A":"CDRL3",
"L95B":"CDRL3",
"L95C":"CDRL3",
"L95D":"CDRL3",
"L95E":"CDRL3",
"L95F":"CDRL3",
"L96":"CDRL3",
"L97":"CDRL3",
"L98":"fwL4",
"L99":"fwL4",
"L100":"fwL4",
"L101":"fwL4",
"L102":"fwL4",
"L103":"fwL4",
"L104":"fwL4",
"L105":"fwL4",
"L106":"fwL4",
"L106A":"fwL4",
"L107":"fwL4",
"L108":"fwL4",
"L109":"fwL4",
"L110":"fwL4",
"L111":"fwL4"}#: abnum format to kabat regions ( framework region 3 divided into two sections a and b)

# Kabat definition for CDRs 
kabat_regions = {}#: Positions of the antibody grouped into region of the structure using abnum format annotations. Chothia definitions for the loops.
for p in positions:
    try:
        kabat_regions[kabat_region_annotations[p]].append(p)
    except KeyError:
        kabat_regions[kabat_region_annotations[p]] = [p]
         
kabat_regions["loops"] = kabat_regions["CDRH1"]+kabat_regions["CDRH2"] +kabat_regions["CDRH3"]+kabat_regions["CDRL1"]+kabat_regions["CDRL2"]+kabat_regions["CDRL3"]
kabat_regions["framework"] = kabat_regions["fwH1"]+kabat_regions["fwH2"]+kabat_regions["afwH3"]+kabat_regions["bfwH3"]+kabat_regions["fwH4"]+kabat_regions["fwL1"]+kabat_regions["fwL2"]+kabat_regions["afwL3"]+kabat_regions["bfwL3"]+kabat_regions["fwL4"] 
# interface positions from 70% non-redundant set, 4.5A any atom contact cutoff and 5% filter of contact frequency.
kabat_regions["interface"] = ["H100", "H100A", "H100B", "H100C", "H100D", "H100E", "H100F", "H100G", "H101", "H102", "H103", "H104", "H105", "H106", "H35", "H37", "H39", "H43", "H44", "H45", "H46", "H47", "H50", "H58", "H59", "H60", "H61", "H62", "H91", "H95", "H96", "H98", "H99", "L1", "L100", "L32", "L34", "L36", "L38", "L41", "L42", "L43", "L44", "L45", "L46", "L49", "L50", "L55", "L85", "L87", "L89", "L91", "L94", "L95", "L95A", "L95B", "L95C", "L96", "L98", "L99"]




# Contact definition for CDRs
contact_region_annotations={"H1":"fwH1",
"H2":"fwH1",
"H3":"fwH1",
"H4":"fwH1",
"H5":"fwH1",
"H6":"fwH1",
"H6A":"fwH1",
"H6B":"fwH1",
"H6C":"fwH1",
"H6D":"fwH1",
"H7":"fwH1",
"H8":"fwH1",
"H9":"fwH1",
"H10":"fwH1",
"H11":"fwH1",
"H12":"fwH1",
"H13":"fwH1",
"H14":"fwH1",
"H15":"fwH1",
"H16":"fwH1",
"H17":"fwH1",
"H18":"fwH1",
"H19":"fwH1",
"H20":"fwH1",
"H21":"fwH1",
"H22":"fwH1",
"H23":"fwH1",
"H24":"fwH1",
"H25":"fwH1",
"H26":"fwH1",
"H27":"fwH1",
"H28":"fwH1",
"H29":"fwH1",
"H30":"CDRH1",
"H31":"CDRH1",
"H31A":"CDRH1",
"H31B":"CDRH1",
"H31C":"CDRH1",
"H31D":"CDRH1",
"H31E":"CDRH1",
"H31F":"CDRH1",
"H31G":"CDRH1",
"H31H":"CDRH1",
"H32":"CDRH1",
"H33":"CDRH1",
"H34":"CDRH1",
"H35":"CDRH1",
"H36":"fwH2",
"H37":"fwH2",
"H38":"fwH2",
"H39":"fwH2",
"H40":"fwH2",
"H41":"fwH2",
"H42":"fwH2",
"H43":"fwH2",
"H44":"fwH2",
"H45":"fwH2",
"H46":"fwH2",
"H47":"CDRH2",
"H48":"CDRH2",
"H49":"CDRH2",
"H50":"CDRH2",
"H51":"CDRH2",
"H52":"CDRH2",
"H52A":"CDRH2",
"H52B":"CDRH2",
"H52C":"CDRH2",
"H53":"CDRH2",
"H54":"CDRH2",
"H55":"CDRH2",
"H56":"CDRH2",
"H57":"CDRH2",
"H58":"CDRH2",
"H59":"afwH3",
"H60":"afwH3",
"H61":"afwH3",
"H62":"afwH3",
"H63":"afwH3",
"H64":"afwH3",
"H65":"afwH3",
"H66":"afwH3",
"H67":"afwH3",
"H68":"afwH3",
"H69":"afwH3",
"H70":"afwH3",
"H71":"afwH3",
"H72":"afwH3",
"H73":"afwH3",
"H74":"afwH3",
"H75":"afwH3",
"H76":"bfwH3",
"H77":"bfwH3",
"H78":"bfwH3",
"H79":"bfwH3",
"H80":"bfwH3",
"H81":"bfwH3",
"H82":"bfwH3",
"H82A":"bfwH3",
"H82B":"bfwH3",
"H82C":"bfwH3",
"H83":"bfwH3",
"H84":"bfwH3",
"H85":"bfwH3",
"H86":"bfwH3",
"H87":"bfwH3",
"H88":"bfwH3",
"H89":"bfwH3",
"H90":"bfwH3",
"H91":"bfwH3",
"H92":"bfwH3",
"H93":"CDRH3",
"H94":"CDRH3",
"H95":"CDRH3",
"H96":"CDRH3",
"H97":"CDRH3",
"H98":"CDRH3",
"H99":"CDRH3",
"H100":"CDRH3",
"H100A":"CDRH3",
"H100B":"CDRH3",
"H100C":"CDRH3",
"H100D":"CDRH3",
"H100E":"CDRH3",
"H100F":"CDRH3",
"H100G":"CDRH3",
"H100H":"CDRH3",
"H100I":"CDRH3",
"H100J":"CDRH3",
"H100K":"CDRH3",
"H100L":"CDRH3",
"H100M":"CDRH3",
"H100N":"CDRH3",
"H100O":"CDRH3",
"H100P":"CDRH3",
"H100Q":"CDRH3",
"H100R":"CDRH3",
"H100S":"CDRH3",
"H100T":"CDRH3",
"H101":"CDRH3",
"H102":"fwH4",
"H103":"fwH4",
"H104":"fwH4",
"H105":"fwH4",
"H106":"fwH4",
"H107":"fwH4",
"H108":"fwH4",
"H109":"fwH4",
"H110":"fwH4",
"H111":"fwH4",
"H112":"fwH4",
"H113":"fwH4",
"L1":"fwL1",
"L2":"fwL1",
"L3":"fwL1",
"L4":"fwL1",
"L5":"fwL1",
"L6":"fwL1",
"L7":"fwL1",
"L8":"fwL1",
"L9":"fwL1",
"L10":"fwL1",
"L11":"fwL1",
"L12":"fwL1",
"L13":"fwL1",
"L14":"fwL1",
"L15":"fwL1",
"L16":"fwL1",
"L17":"fwL1",
"L18":"fwL1",
"L19":"fwL1",
"L20":"fwL1",
"L21":"fwL1",
"L22":"fwL1",
"L23":"fwL1",
"L24":"fwL1",
"L25":"fwL1",
"L26":"fwL1",
"L27":"fwL1",
"L28":"fwL1",
"L29":"fwL1",
"L30":"CDRL1",
"L30A":"CDRL1",
"L30B":"CDRL1",
"L30C":"CDRL1",
"L30D":"CDRL1",
"L30E":"CDRL1",
"L30F":"CDRL1",
"L31":"CDRL1",
"L32":"CDRL1",
"L33":"CDRL1",
"L34":"CDRL1",
"L35":"fwL2",
"L36":"fwL2",
"L37":"fwL2",
"L38":"fwL2",
"L39":"fwL2",
"L40":"fwL2",
"L41":"fwL2",
"L42":"fwL2",
"L43":"fwL2",
"L44":"fwL2",
"L45":"fwL2",
"L46":"CDRL2",
"L47":"CDRL2",
"L48":"CDRL2",
"L49":"CDRL2",
"L50":"CDRL2",
"L51":"CDRL2",
"L52":"CDRL2",
"L53":"CDRL2",
"L54":"CDRL2",
"L54A":"CDRL2",
"L54B":"CDRL2",
"L54C":"CDRL2",
"L54D":"CDRL2",
"L54E":"CDRL2",
"L55":"CDRL2",
"L56":"afwL3",
"L57":"afwL3",
"L58":"afwL3",
"L59":"afwL3",
"L60":"afwL3",
"L61":"afwL3",
"L62":"afwL3",
"L63":"afwL3",
"L64":"afwL3",
"L65":"afwL3",
"L66":"afwL3",
"L66A":"afwL3",
"L66B":"afwL3",
"L66C":"afwL3",
"L66D":"afwL3",
"L66E":"afwL3",
"L66F":"afwL3",
"L66G":"afwL3",
"L66H":"afwL3",
"L67":"afwL3",
"L68":"afwL3",
"L69":"afwL3",
"L70":"afwL3",
"L71":"afwL3",
"L72":"afwL3",
"L73":"bfwL3",
"L74":"bfwL3",
"L75":"bfwL3",
"L76":"bfwL3",
"L77":"bfwL3",
"L78":"bfwL3",
"L79":"bfwL3",
"L80":"bfwL3",
"L81":"bfwL3",
"L82":"bfwL3",
"L83":"bfwL3",
"L84":"bfwL3",
"L85":"bfwL3",
"L86":"bfwL3",
"L87":"bfwL3",
"L88":"bfwL3",
"L89":"bfwL3",
"L90":"bfwL3",
"L91":"bfwL3",
"L92":"bfwL3",
"L93":"CDRL3",
"L94":"CDRL3",
"L95":"CDRL3",
"L95A":"CDRL3",
"L95B":"CDRL3",
"L95C":"CDRL3",
"L95D":"CDRL3",
"L95E":"CDRL3",
"L95F":"CDRL3",
"L96":"CDRL3",
"L97":"fwL4",
"L98":"fwL4",
"L99":"fwL4",
"L100":"fwL4",
"L101":"fwL4",
"L102":"fwL4",
"L103":"fwL4",
"L104":"fwL4",
"L105":"fwL4",
"L106":"fwL4",
"L106A":"fwL4",
"L107":"fwL4",
"L108":"fwL4",
"L109":"fwL4",
"L110":"fwL4",
"L111":"fwL4"}#: abnum format to contact regions ( framework region 3 divided into two sections a and b)

# Contact definition for CDRs 
contact_regions = {}#: Positions of the antibody grouped into region of the structure using abnum format annotations. Chothia definitions for the loops.
for p in positions:
    try:
        contact_regions[contact_region_annotations[p]].append(p)
    except KeyError:
        contact_regions[contact_region_annotations[p]] = [p]
         
contact_regions["loops"] = contact_regions["CDRH1"]+contact_regions["CDRH2"] +contact_regions["CDRH3"]+contact_regions["CDRL1"]+contact_regions["CDRL2"]+contact_regions["CDRL3"]
contact_regions["framework"] = contact_regions["fwH1"]+contact_regions["fwH2"]+contact_regions["afwH3"]+contact_regions["bfwH3"]+contact_regions["fwH4"]+contact_regions["fwL1"]+contact_regions["fwL2"]+contact_regions["afwL3"]+contact_regions["bfwL3"]+contact_regions["fwL4"] 
# interface positions from 70% non-redundant set, 4.5A any atom contact cutoff and 5% filter of contact frequency.
contact_regions["interface"] = ["H100", "H100A", "H100B", "H100C", "H100D", "H100E", "H100F", "H100G", "H101", "H102", "H103", "H104", "H105", "H106", "H35", "H37", "H39", "H43", "H44", "H45", "H46", "H47", "H50", "H58", "H59", "H60", "H61", "H62", "H91", "H95", "H96", "H98", "H99", "L1", "L100", "L32", "L34", "L36", "L38", "L41", "L42", "L43", "L44", "L45", "L46", "L49", "L50", "L55", "L85", "L87", "L89", "L91", "L94", "L95", "L95A", "L95B", "L95C", "L96", "L98", "L99"]



# imgt definition for CDRs
imgt_region_annotations={"H1":"fwH1",
"H2":"fwH1",
"H3":"fwH1",
"H4":"fwH1",
"H5":"fwH1",
"H6":"fwH1",
"H6A":"fwH1",
"H6B":"fwH1",
"H6C":"fwH1",
"H6D":"fwH1",
"H7":"fwH1",
"H8":"fwH1",
"H9":"fwH1",
"H10":"fwH1",
"H11":"fwH1",
"H12":"fwH1",
"H13":"fwH1",
"H14":"fwH1",
"H15":"fwH1",
"H16":"fwH1",
"H17":"fwH1",
"H18":"fwH1",
"H19":"fwH1",
"H20":"fwH1",
"H21":"fwH1",
"H22":"fwH1",
"H23":"fwH1",
"H24":"fwH1",
"H25":"fwH1",
"H26":"CDRH1",
"H27":"CDRH1",
"H28":"CDRH1",
"H29":"CDRH1",
"H30":"CDRH1",
"H31":"CDRH1",
"H31A":"CDRH1",
"H31B":"CDRH1",
"H31C":"CDRH1",
"H31D":"CDRH1",
"H31E":"CDRH1",
"H31F":"CDRH1",
"H31G":"CDRH1",
"H31H":"CDRH1",
"H32":"CDRH1",
"H33":"CDRH1",
"H34":"fwH2",
"H35":"fwH2",
"H36":"fwH2",
"H37":"fwH2",
"H38":"fwH2",
"H39":"fwH2",
"H40":"fwH2",
"H41":"fwH2",
"H42":"fwH2",
"H43":"fwH2",
"H44":"fwH2",
"H45":"fwH2",
"H46":"fwH2",
"H47":"fwH2",
"H48":"fwH2",
"H49":"fwH2",
"H50":"fwH2",
"H51":"CDRH2",
"H52":"CDRH2",
"H52A":"CDRH2",
"H52B":"CDRH2",
"H52C":"CDRH2",
"H53":"CDRH2",
"H54":"CDRH2",
"H55":"CDRH2",
"H56":"CDRH2",
"H57":"CDRH2",
"H58":"afwH3",
"H59":"afwH3",
"H60":"afwH3",
"H61":"afwH3",
"H62":"afwH3",
"H63":"afwH3",
"H64":"afwH3",
"H65":"afwH3",
"H66":"afwH3",
"H67":"afwH3",
"H68":"afwH3",
"H69":"afwH3",
"H70":"afwH3",
"H71":"afwH3",
"H72":"afwH3",
"H73":"afwH3",
"H74":"afwH3",
"H75":"afwH3",
"H76":"bfwH3",
"H77":"bfwH3",
"H78":"bfwH3",
"H79":"bfwH3",
"H80":"bfwH3",
"H81":"bfwH3",
"H82":"bfwH3",
"H82A":"bfwH3",
"H82B":"bfwH3",
"H82C":"bfwH3",
"H83":"bfwH3",
"H84":"bfwH3",
"H85":"bfwH3",
"H86":"bfwH3",
"H87":"bfwH3",
"H88":"bfwH3",
"H89":"bfwH3",
"H90":"bfwH3",
"H91":"bfwH3",
"H92":"bfwH3",
"H93":"CDRH3",
"H94":"CDRH3",
"H95":"CDRH3",
"H96":"CDRH3",
"H97":"CDRH3",
"H98":"CDRH3",
"H99":"CDRH3",
"H100":"CDRH3",
"H100A":"CDRH3",
"H100B":"CDRH3",
"H100C":"CDRH3",
"H100D":"CDRH3",
"H100E":"CDRH3",
"H100F":"CDRH3",
"H100G":"CDRH3",
"H100H":"CDRH3",
"H100I":"CDRH3",
"H100J":"CDRH3",
"H100K":"CDRH3",
"H100L":"CDRH3",
"H100M":"CDRH3",
"H100N":"CDRH3",
"H100O":"CDRH3",
"H100P":"CDRH3",
"H100Q":"CDRH3",
"H100R":"CDRH3",
"H100S":"CDRH3",
"H100T":"CDRH3",
"H101":"CDRH3",
"H102":"CDRH3",
"H103":"fwH4",
"H104":"fwH4",
"H105":"fwH4",
"H106":"fwH4",
"H107":"fwH4",
"H108":"fwH4",
"H109":"fwH4",
"H110":"fwH4",
"H111":"fwH4",
"H112":"fwH4",
"H113":"fwH4",
"L1":"fwL1",
"L2":"fwL1",
"L3":"fwL1",
"L4":"fwL1",
"L5":"fwL1",
"L6":"fwL1",
"L7":"fwL1",
"L8":"fwL1",
"L9":"fwL1",
"L10":"fwL1",
"L11":"fwL1",
"L12":"fwL1",
"L13":"fwL1",
"L14":"fwL1",
"L15":"fwL1",
"L16":"fwL1",
"L17":"fwL1",
"L18":"fwL1",
"L19":"fwL1",
"L20":"fwL1",
"L21":"fwL1",
"L22":"fwL1",
"L23":"fwL1",
"L24":"fwL1",
"L25":"fwL1",
"L26":"fwL1",
"L27":"CDRL1",
"L28":"CDRL1",
"L29":"CDRL1",
"L30":"CDRL1",
"L30A":"CDRL1",
"L30B":"CDRL1",
"L30C":"CDRL1",
"L30D":"CDRL1",
"L30E":"CDRL1",
"L30F":"CDRL1",
"L31":"CDRL1",
"L32":"CDRL1",
"L33":"fwL2",
"L34":"fwL2",
"L35":"fwL2",
"L36":"fwL2",
"L37":"fwL2",
"L38":"fwL2",
"L39":"fwL2",
"L40":"fwL2",
"L41":"fwL2",
"L42":"fwL2",
"L43":"fwL2",
"L44":"fwL2",
"L45":"fwL2",
"L46":"fwL2",
"L47":"fwL2",
"L48":"fwL2",
"L49":"fwL2",
"L50":"CDRL2",
"L51":"CDRL2",
"L52":"CDRL2",
"L53":"afwL3",
"L54":"afwL3",
"L54A":"afwL3",
"L54B":"afwL3",
"L54C":"afwL3",
"L54D":"afwL3",
"L54E":"afwL3",
"L55":"afwL3",
"L56":"afwL3",
"L57":"afwL3",
"L58":"afwL3",
"L59":"afwL3",
"L60":"afwL3",
"L61":"afwL3",
"L62":"afwL3",
"L63":"afwL3",
"L64":"afwL3",
"L65":"afwL3",
"L66":"afwL3",
"L66A":"afwL3",
"L66B":"afwL3",
"L66C":"afwL3",
"L66D":"afwL3",
"L66E":"afwL3",
"L66F":"afwL3",
"L66G":"afwL3",
"L66H":"afwL3",
"L67":"afwL3",
"L68":"afwL3",
"L69":"afwL3",
"L70":"afwL3",
"L71":"afwL3",
"L72":"afwL3",
"L73":"bfwL3",
"L74":"bfwL3",
"L75":"bfwL3",
"L76":"bfwL3",
"L77":"bfwL3",
"L78":"bfwL3",
"L79":"bfwL3",
"L80":"bfwL3",
"L81":"bfwL3",
"L82":"bfwL3",
"L83":"bfwL3",
"L84":"bfwL3",
"L85":"bfwL3",
"L86":"bfwL3",
"L87":"bfwL3",
"L88":"bfwL3",
"L89":"CDRL3",
"L90":"CDRL3",
"L91":"CDRL3",
"L92":"CDRL3",
"L93":"CDRL3",
"L94":"CDRL3",
"L95":"CDRL3",
"L95A":"CDRL3",
"L95B":"CDRL3",
"L95C":"CDRL3",
"L95D":"CDRL3",
"L95E":"CDRL3",
"L95F":"CDRL3",
"L96":"CDRL3",
"L97":"CDRL3",
"L98":"fwL4",
"L99":"fwL4",
"L100":"fwL4",
"L101":"fwL4",
"L102":"fwL4",
"L103":"fwL4",
"L104":"fwL4",
"L105":"fwL4",
"L106":"fwL4",
"L106A":"fwL4",
"L107":"fwL4",
"L108":"fwL4",
"L109":"fwL4",
"L110":"fwL4",
"L111":"fwL4"}#: abnum format to contact regions ( framework region 3 divided into two sections a and b)

# imgt definition for CDRs 
imgt_regions = {}#: Positions of the antibody grouped into region of the structure using abnum format annotations. Chothia definitions for the loops.
for p in positions:
    try:
        imgt_regions[imgt_region_annotations[p]].append(p)
    except KeyError:
        imgt_regions[imgt_region_annotations[p]] = [p]
         
imgt_regions["loops"] = imgt_regions["CDRH1"]+imgt_regions["CDRH2"] +imgt_regions["CDRH3"]+imgt_regions["CDRL1"]+imgt_regions["CDRL2"]+imgt_regions["CDRL3"]
imgt_regions["framework"] = imgt_regions["fwH1"]+imgt_regions["fwH2"]+imgt_regions["afwH3"]+imgt_regions["bfwH3"]+imgt_regions["fwH4"]+imgt_regions["fwL1"]+imgt_regions["fwL2"]+imgt_regions["afwL3"]+imgt_regions["bfwL3"]+contact_regions["fwL4"] 
# interface positions from 70% non-redundant set, 4.5A any atom contact cutoff and 5% filtimgtontact frequency.
imgt_regions["interface"] = ["H100", "H100A", "H100B", "H100C", "H100D", "H100E", "H100F", "H100G", "H101", "H102", "H103", "H104", "H105", "H106", "H35", "H37", "H39", "H43", "H44", "H45", "H46", "H47", "H50", "H58", "H59", "H60", "H61", "H62", "H91", "H95", "H96", "H98", "H99", "L1", "L100", "L32", "L34", "L36", "L38", "L41", "L42", "L43", "L44", "L45", "L46", "L49", "L50", "L55", "L85", "L87", "L89", "L91", "L94", "L95", "L95A", "L95B", "L95C", "L96", "L98", "L99"]



# north definition for CDRs
north_region_annotations={"H1":"fwH1",
"H2":"fwH1",
"H3":"fwH1",
"H4":"fwH1",
"H5":"fwH1",
"H6":"fwH1",
"H6A":"fwH1",
"H6B":"fwH1",
"H6C":"fwH1",
"H6D":"fwH1",
"H7":"fwH1",
"H8":"fwH1",
"H9":"fwH1",
"H10":"fwH1",
"H11":"fwH1",
"H12":"fwH1",
"H13":"fwH1",
"H14":"fwH1",
"H15":"fwH1",
"H16":"fwH1",
"H17":"fwH1",
"H18":"fwH1",
"H19":"fwH1",
"H20":"fwH1",
"H21":"fwH1",
"H22":"fwH1",
"H23":"CDRH1",
"H24":"CDRH1",
"H25":"CDRH1",
"H26":"CDRH1",
"H27":"CDRH1",
"H28":"CDRH1",
"H29":"CDRH1",
"H30":"CDRH1",
"H31":"CDRH1",
"H31A":"CDRH1",
"H31B":"CDRH1",
"H31C":"CDRH1",
"H31D":"CDRH1",
"H31E":"CDRH1",
"H31F":"CDRH1",
"H31G":"CDRH1",
"H31H":"CDRH1",
"H32":"CDRH1",
"H33":"CDRH1",
"H34":"CDRH1",
"H35":"CDRH1",
"H36":"fwH2",
"H37":"fwH2",
"H38":"fwH2",
"H39":"fwH2",
"H40":"fwH2",
"H41":"fwH2",
"H42":"fwH2",
"H43":"fwH2",
"H44":"fwH2",
"H45":"fwH2",
"H46":"fwH2",
"H47":"fwH2",
"H48":"fwH2",
"H49":"fwH2",
"H50":"CDRH2",
"H51":"CDRH2",
"H52":"CDRH2",
"H52A":"CDRH2",
"H52B":"CDRH2",
"H52C":"CDRH2",
"H53":"CDRH2",
"H54":"CDRH2",
"H55":"CDRH2",
"H56":"CDRH2",
"H57":"CDRH2",
"H58":"CDRH2",
"H59":"afwH3",
"H60":"afwH3",
"H61":"afwH3",
"H62":"afwH3",
"H63":"afwH3",
"H64":"afwH3",
"H65":"afwH3",
"H66":"afwH3",
"H67":"afwH3",
"H68":"afwH3",
"H69":"afwH3",
"H70":"afwH3",
"H71":"afwH3",
"H72":"afwH3",
"H73":"afwH3",
"H74":"afwH3",
"H75":"afwH3",
"H76":"bfwH3",
"H77":"bfwH3",
"H78":"bfwH3",
"H79":"bfwH3",
"H80":"bfwH3",
"H81":"bfwH3",
"H82":"bfwH3",
"H82A":"bfwH3",
"H82B":"bfwH3",
"H82C":"bfwH3",
"H83":"bfwH3",
"H84":"bfwH3",
"H85":"bfwH3",
"H86":"bfwH3",
"H87":"bfwH3",
"H88":"bfwH3",
"H89":"bfwH3",
"H90":"bfwH3",
"H91":"bfwH3",
"H92":"bfwH3",
"H93":"CDRH3",
"H94":"CDRH3",
"H95":"CDRH3",
"H96":"CDRH3",
"H97":"CDRH3",
"H98":"CDRH3",
"H99":"CDRH3",
"H100":"CDRH3",
"H100A":"CDRH3",
"H100B":"CDRH3",
"H100C":"CDRH3",
"H100D":"CDRH3",
"H100E":"CDRH3",
"H100F":"CDRH3",
"H100G":"CDRH3",
"H100H":"CDRH3",
"H100I":"CDRH3",
"H100J":"CDRH3",
"H100K":"CDRH3",
"H100L":"CDRH3",
"H100M":"CDRH3",
"H100N":"CDRH3",
"H100O":"CDRH3",
"H100P":"CDRH3",
"H100Q":"CDRH3",
"H100R":"CDRH3",
"H100S":"CDRH3",
"H100T":"CDRH3",
"H101":"CDRH3",
"H102":"CDRH3",
"H103":"fwH4",
"H104":"fwH4",
"H105":"fwH4",
"H106":"fwH4",
"H107":"fwH4",
"H108":"fwH4",
"H109":"fwH4",
"H110":"fwH4",
"H111":"fwH4",
"H112":"fwH4",
"H113":"fwH4",
"L1":"fwL1",
"L2":"fwL1",
"L3":"fwL1",
"L4":"fwL1",
"L5":"fwL1",
"L6":"fwL1",
"L7":"fwL1",
"L8":"fwL1",
"L9":"fwL1",
"L10":"fwL1",
"L11":"fwL1",
"L12":"fwL1",
"L13":"fwL1",
"L14":"fwL1",
"L15":"fwL1",
"L16":"fwL1",
"L17":"fwL1",
"L18":"fwL1",
"L19":"fwL1",
"L20":"fwL1",
"L21":"fwL1",
"L22":"fwL1",
"L23":"fwL1",
"L24":"CDRL1",
"L25":"CDRL1",
"L26":"CDRL1",
"L27":"CDRL1",
"L28":"CDRL1",
"L29":"CDRL1",
"L30":"CDRL1",
"L30A":"CDRL1",
"L30B":"CDRL1",
"L30C":"CDRL1",
"L30D":"CDRL1",
"L30E":"CDRL1",
"L30F":"CDRL1",
"L31":"CDRL1",
"L32":"CDRL1",
"L33":"CDRL1",
"L34":"CDRL1",
"L35":"fwL2",
"L36":"fwL2",
"L37":"fwL2",
"L38":"fwL2",
"L39":"fwL2",
"L40":"fwL2",
"L41":"fwL2",
"L42":"fwL2",
"L43":"fwL2",
"L44":"fwL2",
"L45":"fwL2",
"L46":"fwL2",
"L47":"fwL2",
"L48":"fwL2",
"L49":"CDRL2",
"L50":"CDRL2",
"L51":"CDRL2",
"L52":"CDRL2",
"L53":"CDRL2",
"L54":"CDRL2",
"L54A":"CDRL2",
"L54B":"CDRL2",
"L54C":"CDRL2",
"L54D":"CDRL2",
"L54E":"CDRL2",
"L55":"CDRL2",
"L56":"CDRL2",
"L57":"afwL3",
"L58":"afwL3",
"L59":"afwL3",
"L60":"afwL3",
"L61":"afwL3",
"L62":"afwL3",
"L63":"afwL3",
"L64":"afwL3",
"L65":"afwL3",
"L66":"afwL3",
"L66A":"afwL3",
"L66B":"afwL3",
"L66C":"afwL3",
"L66D":"afwL3",
"L66E":"afwL3",
"L66F":"afwL3",
"L66G":"afwL3",
"L66H":"afwL3",
"L67":"afwL3",
"L68":"afwL3",
"L69":"afwL3",
"L70":"afwL3",
"L71":"afwL3",
"L72":"afwL3",
"L73":"bfwL3",
"L74":"bfwL3",
"L75":"bfwL3",
"L76":"bfwL3",
"L77":"bfwL3",
"L78":"bfwL3",
"L79":"bfwL3",
"L80":"bfwL3",
"L81":"bfwL3",
"L82":"bfwL3",
"L83":"bfwL3",
"L84":"bfwL3",
"L85":"bfwL3",
"L86":"bfwL3",
"L87":"bfwL3",
"L88":"bfwL3",
"L89":"CDRL3",
"L90":"CDRL3",
"L91":"CDRL3",
"L92":"CDRL3",
"L93":"CDRL3",
"L94":"CDRL3",
"L95":"CDRL3",
"L95A":"CDRL3",
"L95B":"CDRL3",
"L95C":"CDRL3",
"L95D":"CDRL3",
"L95E":"CDRL3",
"L95F":"CDRL3",
"L96":"CDRL3",
"L97":"CDRL3",
"L98":"fwL4",
"L99":"fwL4",
"L100":"fwL4",
"L101":"fwL4",
"L102":"fwL4",
"L103":"fwL4",
"L104":"fwL4",
"L105":"fwL4",
"L106":"fwL4",
"L106A":"fwL4",
"L107":"fwL4",
"L108":"fwL4",
"L109":"fwL4",
"L110":"fwL4",
"L111":"fwL4"}#: abnum format to contact regions ( framework region 3 divided into two sections a and b)

# north definition for CDRs 
north_regions = {}#: Positions of the antibody grouped into region of the structure using abnum format annotations. Chothia definitions for the loops.
for p in positions:
    try:
        north_regions[north_region_annotations[p]].append(p)
    except KeyError:
        north_regions[north_region_annotations[p]] = [p]
         
north_regions["loops"] = north_regions["CDRH1"]+north_regions["CDRH2"] +north_regions["CDRH3"]+north_regions["CDRL1"]+north_regions["CDRL2"]+north_regions["CDRL3"]
north_regions["framework"] = north_regions["fwH1"]+north_regions["fwH2"]+north_regions["afwH3"]+north_regions["bfwH3"]+north_regions["fwH4"]+north_regions["fwL1"]+north_regions["fwL2"]+north_regions["afwL3"]+north_regions["bfwL3"]+contact_regions["fwL4"] 
# interface positions from 70% non-redundant set, 4.5A any atom contact cutoff and 5% filtimgtontact frequency.
north_regions["interface"] = ["H100", "H100A", "H100B", "H100C", "H100D", "H100E", "H100F", "H100G", "H101", "H102", "H103", "H104", "H105", "H106", "H35", "H37", "H39", "H43", "H44", "H45", "H46", "H47", "H50", "H58", "H59", "H60", "H61", "H62", "H91", "H95", "H96", "H98", "H99", "L1", "L100", "L32", "L34", "L36", "L38", "L41", "L42", "L43", "L44", "L45", "L46", "L49", "L50", "L55", "L85", "L87", "L89", "L91", "L94", "L95", "L95A", "L95B", "L95C", "L96", "L98", "L99"]



sci_to_common_names={"camelus dromedarius" : "arabian camel",
"dermatophagoides farinae" : "american house dust mite",
"phasianus colchicus" : "ring-necked pheasant",
"oryctolagus cuniculus" : "rabbit",
"lama pacos" : "alpaca",
"saccharomyces cerevisiae" : "baker's yeast",
"numida meleagris" : "helmeted guineafowl",
"ovis aries" : "sheep",
"blattella germanica" : "german cockroach",
"rattus norvegicus" : "norway rat",
"ictalurus punctatus" : "channel catfish",
"phoca vitulina" : "harbor seal",
"rattus rattus" : "black rat",
"sus scrofa" : "pig",
"centruroides noxius hoffmann" : "mexican scorpion",
"betula pendula" : "european white birch",
"androctonus australis hector" : "fat-tailed scorpion",
"llama glama" : "llama",
"camelidae" : "camelid",
"macaca mulatta" : "rhesus monkey",
"phleum pratense" : "common timothy",
"influenza a virus" : "influenza a virus",
"homo sapiens" : "human",
"colinus virginianus" : "northern bobwhite",
"dermatophagoides pteronyssinus" : "european house dust mite",
"equus caballus" : "horse",
"caenorhabditis elegans" : "nematode",
"meleagris gallopavo" : "turkey",
"bos taurus" : "cattle",
"mus musculus" : "house mouse",
"gallus gallus" : "chicken",
"apis mellifera" : "honey bee",
"lama glama" : "llama",
"vicugna pacos" : "alpaca",
"ricinus communis" : "castor bean",
"aequorea victoria" : "jellyfish",
"cricetulus migratorius" : "armenian hamster",
"pan troglodytes" : "chimpanzee"}

def tuple_interpret(x):
    """
    Interpretation for numbering from abnum into tuple numbering. 
    
    @param x: A string for the identifier of a numbered position. e.g "H100A". 
    @type x: C{str}
    
    @return: A tuple of the chain tupe followed by a tuple of residue id and insertion code. eg. ( H, (100, "A") )
    @rtype: C{tuple}
    
    """
    assert x[0] == "H" or x[0] == "L", x
    try:
        return x[0],( int(x[1:]), ' ')
    except ValueError:
        return x[0],( int(x[1:-1]), x[-1] )


chothia_regions_tuples = {"H":{},"L":{}}#: Positions of the antibody grouped into regions of the structure using Bio.PDB format annotations. Chothia definitions for the loops.
for p in positions:
    c,t = tuple_interpret(p)
    try:
        chothia_regions_tuples[c][chothia_region_annotations[p]].append(t)
    except KeyError:
        chothia_regions_tuples[c][chothia_region_annotations[p]] = [t]
    if p in regions["interface"]:
        try:
            chothia_regions_tuples[c]["interface"].append(t)
        except KeyError:
            chothia_regions_tuples[c]["interface"]=[t]

chothia_regions_tuples["H"]["loops"] = chothia_regions_tuples["H"]["CDRH1"]+chothia_regions_tuples["H"]["CDRH2"] +chothia_regions_tuples["H"]["CDRH3"]
chothia_regions_tuples["L"]["loops"] = chothia_regions_tuples["L"]["CDRL1"]+chothia_regions_tuples["L"]["CDRL2"] +chothia_regions_tuples["L"]["CDRL3"]
chothia_regions_tuples["H"]["framework"] = chothia_regions_tuples["H"]["fwH1"]+chothia_regions_tuples["H"]["fwH2"]+chothia_regions_tuples["H"]["afwH3"]+chothia_regions_tuples["H"]["bfwH3"]+chothia_regions_tuples["H"]["fwH4"]
chothia_regions_tuples["L"]["framework"] = chothia_regions_tuples["L"]["fwL1"]+chothia_regions_tuples["L"]["fwL2"]+chothia_regions_tuples["L"]["afwL3"]+chothia_regions_tuples["L"]["bfwL3"]+chothia_regions_tuples["L"]["fwL4"] 

regions_tuples = chothia_regions_tuples



kabat_regions_tuples = {"H":{},"L":{}}#: Positions of the antibody grouped into regions of the structure using Bio.PDB format annotations. kabat definitions for the loops.
for p in positions:
    c,t = tuple_interpret(p)
    try:
        kabat_regions_tuples[c][kabat_region_annotations[p]].append(t)
    except KeyError:
        kabat_regions_tuples[c][kabat_region_annotations[p]] = [t]
    if p in regions["interface"]:
        try:
            kabat_regions_tuples[c]["interface"].append(t)
        except KeyError:
            kabat_regions_tuples[c]["interface"]=[t]

kabat_regions_tuples["H"]["loops"] = kabat_regions_tuples["H"]["CDRH1"]+kabat_regions_tuples["H"]["CDRH2"] +kabat_regions_tuples["H"]["CDRH3"]
kabat_regions_tuples["L"]["loops"] = kabat_regions_tuples["L"]["CDRL1"]+kabat_regions_tuples["L"]["CDRL2"] +kabat_regions_tuples["L"]["CDRL3"]
kabat_regions_tuples["H"]["framework"] = kabat_regions_tuples["H"]["fwH1"]+kabat_regions_tuples["H"]["fwH2"]+kabat_regions_tuples["H"]["afwH3"]+kabat_regions_tuples["H"]["bfwH3"]+kabat_regions_tuples["H"]["fwH4"]
kabat_regions_tuples["L"]["framework"] = kabat_regions_tuples["L"]["fwL1"]+kabat_regions_tuples["L"]["fwL2"]+kabat_regions_tuples["L"]["afwL3"]+kabat_regions_tuples["L"]["bfwL3"]+kabat_regions_tuples["L"]["fwL4"] 


contact_regions_tuples = {"H":{},"L":{}}#: Positions of the antibody grouped into regions of the structure using Bio.PDB format annotations. contact definitions for the loops.
for p in positions:
    c,t = tuple_interpret(p)
    try:
        contact_regions_tuples[c][contact_region_annotations[p]].append(t)
    except KeyError:
        contact_regions_tuples[c][contact_region_annotations[p]] = [t]
    if p in regions["interface"]:
        try:
            contact_regions_tuples[c]["interface"].append(t)
        except KeyError:
            contact_regions_tuples[c]["interface"]=[t]

contact_regions_tuples["H"]["loops"] = contact_regions_tuples["H"]["CDRH1"]+contact_regions_tuples["H"]["CDRH2"] +contact_regions_tuples["H"]["CDRH3"]
contact_regions_tuples["L"]["loops"] = contact_regions_tuples["L"]["CDRL1"]+contact_regions_tuples["L"]["CDRL2"] +contact_regions_tuples["L"]["CDRL3"]
contact_regions_tuples["H"]["framework"] = contact_regions_tuples["H"]["fwH1"]+contact_regions_tuples["H"]["fwH2"]+contact_regions_tuples["H"]["afwH3"]+contact_regions_tuples["H"]["bfwH3"]+contact_regions_tuples["H"]["fwH4"]
contact_regions_tuples["L"]["framework"] = contact_regions_tuples["L"]["fwL1"]+contact_regions_tuples["L"]["fwL2"]+contact_regions_tuples["L"]["afwL3"]+contact_regions_tuples["L"]["bfwL3"]+contact_regions_tuples["L"]["fwL4"] 

imgt_regions_tuples = {"H":{},"L":{}}#: Positions of the antibody grouped into regions of the structure using Bio.PDB format annotations. imgt definitions for the loops.
for p in positions:
    c,t = tuple_interpret(p)
    try:
        imgt_regions_tuples[c][imgt_region_annotations[p]].append(t)
    except KeyError:
        imgt_regions_tuples[c][imgt_region_annotations[p]] = [t]
    if p in regions["interface"]:
        try:
            imgt_regions_tuples[c]["interface"].append(t)
        except KeyError:
            imgt_regions_tuples[c]["interface"]=[t]

imgt_regions_tuples["H"]["loops"] = imgt_regions_tuples["H"]["CDRH1"]+imgt_regions_tuples["H"]["CDRH2"] +imgt_regions_tuples["H"]["CDRH3"]
imgt_regions_tuples["L"]["loops"] = imgt_regions_tuples["L"]["CDRL1"]+imgt_regions_tuples["L"]["CDRL2"] +imgt_regions_tuples["L"]["CDRL3"]
imgt_regions_tuples["H"]["framework"] = imgt_regions_tuples["H"]["fwH1"]+imgt_regions_tuples["H"]["fwH2"]+imgt_regions_tuples["H"]["afwH3"]+imgt_regions_tuples["H"]["bfwH3"]+imgt_regions_tuples["H"]["fwH4"]
imgt_regions_tuples["L"]["framework"] = imgt_regions_tuples["L"]["fwL1"]+imgt_regions_tuples["L"]["fwL2"]+imgt_regions_tuples["L"]["afwL3"]+imgt_regions_tuples["L"]["bfwL3"]+imgt_regions_tuples["L"]["fwL4"] 

north_regions_tuples = {"H":{},"L":{}}#: Positions of the antibody grouped into regions of the structure using Bio.PDB format annotations. north definitions for the loops.
for p in positions:
    c,t = tuple_interpret(p)
    try:
        north_regions_tuples[c][north_region_annotations[p]].append(t)
    except KeyError:
        north_regions_tuples[c][north_region_annotations[p]] = [t]
    if p in regions["interface"]:
        try:
            north_regions_tuples[c]["interface"].append(t)
        except KeyError:
            north_regions_tuples[c]["interface"]=[t]

north_regions_tuples["H"]["loops"] = north_regions_tuples["H"]["CDRH1"]+north_regions_tuples["H"]["CDRH2"] +north_regions_tuples["H"]["CDRH3"]
north_regions_tuples["L"]["loops"] = north_regions_tuples["L"]["CDRL1"]+north_regions_tuples["L"]["CDRL2"] +north_regions_tuples["L"]["CDRL3"]
north_regions_tuples["H"]["framework"] = north_regions_tuples["H"]["fwH1"]+north_regions_tuples["H"]["fwH2"]+north_regions_tuples["H"]["afwH3"]+north_regions_tuples["H"]["bfwH3"]+north_regions_tuples["H"]["fwH4"]
north_regions_tuples["L"]["framework"] = north_regions_tuples["L"]["fwL1"]+north_regions_tuples["L"]["fwL2"]+north_regions_tuples["L"]["afwL3"]+north_regions_tuples["L"]["bfwL3"]+north_regions_tuples["L"]["fwL4"] 

def is_interface_region(entity):
    if hasattr(entity, "region"):
        region = entity.region
        if region == "?": return False
        c = region[-2]
        ident = tuple(entity.id[1:])
        if ident in regions_tuples[c]["interface"]:
            return True
        else:
            return False
    else:
        raise Exception("Unknown region")
        
        
def get_position_region(position,definition="chothia"):
    """
    Get the region annotation for a certain chothia position.
    
    @param position: An abnum position e.g. "H100A"
    @type position: C{str}
    
    @return: The region of the antibody that the position is in. e.g. "CDRH3"
    @rtype: C{str}
    """
    
    definition=definition.lower()
    assert definition in ["chothia","kabat","contact","imgt","north"], "Unrecognised CDR definition: %s"%definition

    return {"chothia":chothia_region_annotations, "kabat":kabat_region_annotations,"contact":contact_region_annotations, 'imgt':imgt_region_annotations, 'north':north_region_annotations}[definition][position] 

def find_identity( seq1 ,seq2, positions=[]):
    """
    Find the matched sequence identity between two aligned sequences.
    
    @param seq1: Dictionary with key as the position and value as the single letter amino acid code. or an aligned list or string
    @param seq2: Dictionary with key as the position and value as the single letter amino acid code. or an aligned list or string
    @param positions: List of positions over which to calculate sequence identity. 
    """

    if isinstance(seq1,dict) and isinstance(seq2,dict):  
        if not positions:
            positions = list(seq1.keys())
    else:
        if not positions:
            positions = list(range(len(seq1)))
        
    n = 0 #  number 
    m = 0 #  match 
    # matched identity
    for p in positions:
        try:
            if seq1[p] == "-":
                continue
            if seq2[p] == "-":
                continue
        except KeyError:
            continue
        
        if seq1[p] == seq2[p]:
            m +=1
        n+=1
    try:
        return float(m)/n
    except ZeroDivisionError:
        print("Warning no common positions given")
        return 0

def uniq(seq, idfun=None):
    """
    A function to uniquify a sequence.
    With thanks to http://www.peterbe.com

    @param seq: A sequence to uniquify
    @param idfun: An optional function to use as a key. Like the "key" kwarg in C{sorted}. 
    
    @return: The sequence.
    """
    # order preserving
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        if marker in seen: continue
        seen[marker] = 1
        result.append(item)
    return result


def maxind(l, MIN=False, value=None): 
    """
    A recursive function to find the indices of the maximum (or minimum) values of a list l.
    
    @param l: A sequence to find the indices of the maximum values.
    @param MIN: Flag to find the min instead of the max
    @param value: An optional argument to find indices of a particular value 
    
    @return: A list of the indices of the list which have the maximum (or minimum or specific) value. 
    @rtype: C{list}
             
    """
    if value:
        m = value
    elif MIN:
        m = min(l)
    else:    
        m = max(l)
    def rec(i):
        try:    
            i.append( l.index( m, i[-1]+1 ) )
        except ValueError:
            return i
        rec(i)
        return i[1:]
    return rec([-1])

def get_coords(fv, regions, atom_types):
    """
    Return a dictionary of coordinates from the structure which correspond to the part you are interested in.  
    regions is a list of keyword value tuples which specify the region required.
    atom_types specify which atom types you want default Calpha
    fv is a structure object.
    """
    picked_coords = {}
    coord_names = []        
    for atom in fv.get_atoms():
        # Query the atoms in the structure.
        if atom.id in atom_types or "all" in atom_types:
            for keyword, value in regions:
                if keyword == "chain_type" and atom.chain_type == value:
                    picked_coords[ (atom.chain_type, atom.parent.id[1], atom.parent.id[2], atom.id ) ] = atom.coord
                    coord_names.append((atom.chain_type, atom.parent.id[1], atom.parent.id[2], atom.id ) )
                    break
                elif keyword == "region" and atom.region == value:
                    picked_coords[ (atom.chain_type, atom.parent.id[1], atom.parent.id[2], atom.id ) ] = atom.coord
                    coord_names.append((atom.chain_type, atom.parent.id[1], atom.parent.id[2], atom.id ) )
                    break
                elif keyword == "residue" and (atom.chain_type,atom.parent.id[1:]) == tuple_interpret(value):
                    picked_coords[ (atom.chain_type, atom.parent.id[1], atom.parent.id[2], atom.id ) ] = atom.coord
                    coord_names.append((atom.chain_type, atom.parent.id[1], atom.parent.id[2], atom.id ) )
                    break
    assert len(picked_coords) > 0, "No coordinates found that satisfy the selection: %s"%repr(regions)
    return picked_coords, coord_names  


def which(name, flags=os.X_OK):
    """
    Search PATH for executable files with the given name.
   
    On newer versions of MS-Windows, the PATHEXT environment variable will be
    set to the list of file extensions for files considered executable. This
    will normally include things like ".EXE". This fuction will also find files
    with the given name ending with any of these extensions.

    On MS-Windows the only flag that has any meaning is os.F_OK. Any other
    flags will be ignored.
   
    @type name: C{str}
    @param name: The name for which to search.
   
    @type flags: C{int}
    @param flags: Arguments to L{os.access}.
   
    @rtype: C{list}
    @param: A list of the unique full paths to files found, in the
    order in which they were found.
    """
    result = []
    exts = [_f for _f in os.environ.get('PATHEXT', '').split(os.pathsep) if _f]
    path = os.environ.get('PATH', None)
    if path is None:
        return []
    for p in os.environ.get('PATH', '').split(os.pathsep):
        p = os.path.join(p, name)
        if os.access(p, flags):
            result.append(p)
        for e in exts:
            pext = p + e
            if os.access(pext, flags):
                result.append(pext)
    return uniq(result)
