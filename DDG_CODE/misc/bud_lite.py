#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys

from klab.fs.fsio import read_file, write_file, get_file_lines, write_temp_file

sys.path.insert(0, "../..")
from ddg.ddglib.ppi_api import get_interface as get_ppi_interface


ppi_api = get_ppi_interface(read_file('pw'))
ppi_api.DDG_db._get_connection()


def generate_bud_lite():
    ids = []
    lite_predictions = ppi_api.DDG_db.execute_select('SELECT ID, UserPPDataSetExperimentID FROM PredictionPPI WHERE PredictionSet="ddg_monomer_16-zemu-betajuly15" AND (maxvmem IS NOT null AND maxvmem <= 3.5 AND ddGTime <= 3500)')

    con = ppi_api.DDG_db.connection
    with con:
        cur = con.cursor()
        for p in lite_predictions:
            ids.append(p['ID'])
            try:
                qry = 'INSERT INTO UserPPDataSetExperimentTag (UserPPDataSetExperimentID, Tag) VALUES({0},"{1}")'.format(p['UserPPDataSetExperimentID'], 'ZEMu lite')
                cur.execute(qry)
            except:
                pass
        cur.execute('SELECT COUNT(Tag) AS NumInserted FROM UserPPDataSetExperimentTag WHERE Tag="ZEMu lite"')
        print('There are {0} cases of bud lite in the fridge #trademarknotrademark.'.format(cur.fetchone()['NumInserted']))
    return ids


def generate_smirnoff_black():
    ids = []
    heavy_fuel_predictions = ppi_api.DDG_db.execute_select('SELECT ID, UserPPDataSetExperimentID FROM PredictionPPI WHERE PredictionSet="ddg_monomer_16-zemu-betajuly15" AND NOT(maxvmem IS NOT null AND maxvmem <= 3.5 AND ddGTime <= 3500)')

    con = ppi_api.DDG_db.connection
    with con:
        cur = con.cursor()
        for p in heavy_fuel_predictions:
            ids.append(p['ID'])
            try:
                qry = 'INSERT INTO UserPPDataSetExperimentTag (UserPPDataSetExperimentID, Tag) VALUES({0},"{1}")'.format(p['UserPPDataSetExperimentID'], 'ZEMu heavy')
                cur.execute(qry)
            except:
                pass
        cur.execute('SELECT COUNT(Tag) AS NumInserted FROM UserPPDataSetExperimentTag WHERE Tag="ZEMu heavy"')
        print('There are {0} cases of heavy, heavy fuel in the fridge #straitsnodirestraits.'.format(cur.fetchone()['NumInserted']))
    return ids


lite_ids = generate_bud_lite()
heavy_souls = generate_smirnoff_black()


all_predictions = set([r['ID'] for r in ppi_api.DDG_db.execute_select('SELECT ID FROM PredictionPPI WHERE PredictionSet="ddg_monomer_16-zemu-betajuly15"')])
assert(len(lite_ids) == len(set(lite_ids)))
assert(len(heavy_souls) == len(set(heavy_souls)))
lite_ids = set(lite_ids)
heavy_souls = set(heavy_souls)
assert(len(lite_ids.intersection(heavy_souls)) == 0)
assert(len(lite_ids.difference(all_predictions)) == 0)
assert(len(heavy_souls.difference(all_predictions)) == 0)
missing_predictions = all_predictions.difference(lite_ids.union(heavy_souls))
assert(len(missing_predictions) == 0)
