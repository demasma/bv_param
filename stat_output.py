#!/usr/bin/python

import pg
import numpy
import string
import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-i', '--input_metal', dest='metal_type', \
                  default = None, help='Provide PDB ID code')
parser.add_option('-t', '--threshold', dest='threshold', \
                  default = None, help='threshold')
(options, args) = parser.parse_args()


atom_type = options.metal_type
threshold = options.threshold
table_1cs = "{0}_param_r0".format(atom_type)
table_4a = "{0}_4a_param_r0".format(atom_type)
BASEDIR = '/home/matt_d/Documents/Lab/Projects/bondvalence_parameters/'
DIRECTORY = BASEDIR + 'figures/'

# establish database connection
bvp = pg.connect('bvparam2012', 'localhost')

# input list, identifier
# output histogram in gnuplot readable format (*.dat file)
def hist_output(inlist, filename, in_bin):
    '''INPUTS AND OUTPUTS
    inlist
    filename'''
    global DIRECTORY
    if not filename.endswith('.dat'):
        filename+='.dat'
    if not os.path.exists(DIRECTORY):
        os.makedirs(DIRECTORY)
    FIGDIR = DIRECTORY + atom_type
    if not os.path.exists(FIGDIR):
        os.makedirs(FIGDIR)
    filename = FIGDIR + os.sep + atom_type + '-' + filename
    outfile = open(filename, 'w+')
    #tmphist = numpy.histogram(inlist,bins=10)
    tmphist = numpy.histogram(inlist, bins = in_bin)
    MAX_HEIGHT = numpy.max(tmphist[0])
    tmpbin = tmphist[1][0:-1]
    bin_width = tmpbin[1] - tmpbin[0]
    bin_center = tmpbin + bin_width 
    for j in range(len(tmphist[0])):
        # uncomment to output normalize bins
        # write_string = "{0}\t{1}\n".format(str(tmphist[1][j]),\
        #                str(float(tmphist[0][i])/float(len(inlist))))
        write_string = "{0}\t{1}\t{2}\n".format(str(bin_center[j]),\
                       str(float(tmphist[0][j])), str(bin_width))
        outfile.write(write_string)
    outfile.close()
    # get the histogram and bins and height of tallest bin
#    filename, bins, pk_cnt = hist_output(mixr0, "{0}_{1}".format(desc,\
#                                             str(anion_num)), in_bin)
#    stats_dict['pk_cnt'] = int(pk_cnt)    
    return [filename, tmphist[1], MAX_HEIGHT]

def get_atom_name(atomic_number):
    query = 'select elem from atom_number where atomnum = \
            {0}'.format(str(atomic_number))
    query = bvp.query(query)
    an_sym = query.getresult()
    an_sym = an_sym[0][0].replace(' ','')
    return an_sym

def gnuplotter(files):
    global atom_type
    global directory
    filename = directory + "{0}_histograms.gp".format(atom_type)
    outfile = open(filename, 'w+')
    plotting_str = "\
                    set terminal png giant truecolor interlace enhanced \
                        size 1920,1080 \n\
                    set encoding iso_8859_1 \n\
                    set key outside \n\
                    set grid \n\
                    set log y \n"
    outfile.write(plotting_str)
#    plot three histograms on one along with the and mean values
    for k in files.keys():
        an_sym = get_atom_number(k)
        outstr = "set output '{3}/{2}_{0}.png'\n\
                  unset arrow\n\
                  unset label\n\
                  set title 'breakdown {2} and {1}'\n\
                  ".format(str(k), an_sym, atom_type,directory)
        outstr2 = 'plot '
        outstr3 = ''
        lc_cnt = 10
        for n in range(len(files[k][0])):
            outstr2 += '''\
                       "{0}" u 1:2:3 w boxes fs solid 0.25 border'\
                       t '{1}' lc {2}\
                       '''.format(files[k][n][0], files[k][n][0].\
                                  replace('.dat', '').split('/')[-1], \
                                  str(lc_cnt))
            outstr3 += "set style arrow 5 nohead front lw 2 lc rgb 'black'\n"\
                           .format(str(lc_cnt))
            outstr3 += "set arrow from {0:.2f}, graph 0 to {0:.2f},graph 1 as 5\n"\
                           .format(files[k][n][1]['median'], \
                                   `float(files[k][n][2]['pk_cnt'])/2.0`)
            outstr3 += "set label 'median' at {0:.2f},{1} right front\
                        textcolor rgb 'red'\n"\
                           .format(files[k][n][1]['median'], \
                                   `float(files[k][n][1]['pk_cnt'])/2.0`)
            outstr3 += "set style arrow 6 nohead empty front lw 2 lc {0} lt 5\n"\
                           .format(str(lc_cnt))
            #outstr3 += "set arrow from {0:.2f},0 to {0:.2f},{1} as 6\n"\
            outstr3 += "set arrow from {0:.2f}, graph 0 to {0:.2f}, graph 0.5\
                            as 6\n"\
                           .format(files[k][n][1]['mean'], \
                                   `float(files[k][n][1]['pk_cnt'])/1.0`)
            if files[k][n][0] is not files[k][-1][0]:
                outstr2 += ", "
            else:
                outstr2 += "\n"
            lc_cnt += 1
        outfile.write(outstr + outstr3 + outstr2)
    outfile.close()

def select_query(info, table, modifiers):
    query_base = "SELECT {0} FROM {1}".format(info, table)
    if type(modifiers) is list:
        i = 0
        for arg in modifiers:
            if i == 0:
                query_base += " {0}".format(arg)
            else:
                query_base += " AND {0}".format(arg)
            i += 1
    elif type(modifiers) is str:
        query_base += " {0}".format(modifiers)
    return query_base

def stat_output(anion_num, desc, percentage, in_bin, table, stat_tab):
    '''anion_num should be a tuple returned from database query
    desc is the output file prefix
    percentage is a string which is apppended to the query selecting upper or 
               lower 50 percentiles'''
    global atom_type
    if atom_type == "iron2":
        cut_dict = {6: ">1.49", 7: ">1.66", 8: ">1.66", 9: ">1.59", 15:">2.06", 16:">2.02", \
                    17:">1.99", 34:">2.16", 35:">2.14", 53:">2.34"}
    elif atom_type == "iron3":
        cut_dict = {6: "<1.64", 7: "<1.81", 8: "<1.80", 9: "<1.74", 15: "<2.21 and valence>2.12", 16: "<2.17 and valence>2.11", \
                    17: "<2.14", 34: "<2.31", 35: "<2.29", 53: "<2.48"}

    # query the database for the proper cation-anion connections
    # and select the appropriate percentage
    query = select_query("valence", table, ['WHERE atomnum = {0}'\
                         .format(str(anion_num)), percentage, \
                "valence{0}".format(str(cut_dict[anion_num]))])
    print query
    mixr0tab = bvp.query(query)
    mixr0 = mixr0tab.getresult()

    # get number of homoleptic sites
    homo_query = select_query("count(*)", table, ["WHERE atomnum = {0}"\
                              .format(str(anion_num)), "percentage = 1", \
                              "valence{0}".format(cut_dict[anion_num])])
    print homo_query                            
    homo_cnt_tab = bvp.query(homo_query)
    homo_cnt = homo_cnt_tab.getresult()[0][0]
    print anion_num, homo_cnt

    # get minimum valence
    min_query = select_query("min(valence)", table, ["WHERE atomnum = {0}"\
                .format(str(anion_num)), percentage])
    min_val_tab = bvp.query(min_query)
    min_val = min_val_tab.getresult()[0][0]

    # get maximum valence
    max_query = select_query("max(valence)", table, ["WHERE atomnum = {0}"\
                .format(str(anion_num)), percentage])
    max_val_tab = bvp.query(max_query)
    max_val = max_val_tab.getresult()[0][0]

    # get average valence
    avg_query = select_query("avg(valence)", table, ["WHERE atomnum = {0}"\
                .format(str(anion_num)), percentage, \
                "valence{0}".format(cut_dict[anion_num])])
    avg_val_tab = bvp.query(avg_query)
    avg_val = avg_val_tab.getresult()[0][0]

    # get stddev of valence
    stddev_query = select_query("stddev(valence)", table, ["WHERE atomnum = \
                                {0}".format(str(anion_num)), percentage, \
                "valence{0}".format(cut_dict[anion_num])])
    stddev_val_tab = bvp.query(stddev_query)
    stddev_val = stddev_val_tab.getresult()[0][0]

    # convert the list of single-valued tuples into a 'vector'
    mixr0 = [item for sublist in mixr0 for item in sublist]
    # get number of negative R_0 values
    query += " AND valence < 0"
    neg_cnt_tab = bvp.query(query)
    neg_cnt = neg_cnt_tab.getresult()

    # compose the dictionary of statistics and prepare for insertion
    # into table

    an_name = get_atom_name(anion_num)
    stats_dict = {'an_num': anion_num, 'total': len(mixr0), \
                  'max': max_val, 'mean': avg_val, \
                  'median': numpy.median(mixr0), 'stddev': numpy.std(mixr0), \
                  'num_neg': len(neg_cnt), 'min': min_val, \
                  'comment': desc, 'num_homo': homo_cnt, 'atom_name': an_name}

    filename = "filename"
    bins = "bins"
    for j in stats_dict:
        stats_dict[j] = str(stats_dict[j])
    bvpdb.insert('public.{0}'.format(stat_tab), stats_dict)
    return [filename, bins, stats_dict]

def stat_table_setup(in_tab_name, percentage, coord_sph):
    # get unique list of atoms
    anion_type_query = select_query('distinct(atomnum)', in_tab_name, \
                      'WHERE percentage >= ' + threshold)
    anion_type_res = bvp.query(anion_type_query)
    anion_list = anion_type_res.getresult()
    anion_list = tup_list_flatten(anion_list)

    out_tab_name = "{0}_{1}_test_stats".format(atom_type, coord_sph)
    bvp.query('DROP TABLE IF EXISTS public.{0}'.format(out_tab_name))
    query = 'CREATE TABLE public.{0} (\
             an_num smallint, total integer, mean double precision, \
             stddev double precision, median double precision, \
             max double precision, min double precision, num_neg integer, \
             num_homo integer, comment text, atom_name text)'\
             .format(out_tab_name)
    bvp.query(query)

    return anion_list, out_tab_name

def tup_list_flatten(tup_list):
    flat_list = [item for sublist in tup_list for item in sublist]
    return flat_list

def stat_printer(stat_dict, key_to_print):
    col_wid = 10
    matrix_output = []
    matrix_output.append(key_to_print)
    #print '\033[1m' + '\033[91m',
    #for i in key_to_print:
    #    print i.rjust(col_wid),
    #print '\033[0m'
    stat_dict_keys = stat_dict.keys()
    for i in sorted(stat_dict_keys):
        for k in range(len(stat_dict[i])):
            tmplist =[]
            #if k == 0:
            #    print '\033[33m',
            #else:
            #    print '\033[34m',
            for j in out_order:
                tmplist.append(stat_dict[i][k][1][j])
            #    #if type(stat_dict[i][k][1][j]) is float:
            #    #    print '{0: .3f}'.format(stat_dict[i][k][1][j]).rjust(col_wid),
            #    #else:
            #    #    print `stat_dict[i][k][1][j]`.rjust(col_wid),
            #print
            matrix_output.append(tmplist)
    #print '\033[0m'
    matrix_output = numpy.matrix(matrix_output)
    trans_matrix = numpy.transpose(matrix_output)
    for row in range(trans_matrix.shape[0]):
        for ele in range(trans_matrix.shape[1]):
            print trans_matrix[row, ele].rjust(11),
        print
if __name__ == "__main__":
    an_list_1cs, stat_tab_1cs = stat_table_setup(table_1cs, threshold, "1cs")
    an_list_4a, stat_tab_4a = stat_table_setup(table_4a, threshold, "4a")

    # create statistics table for given element
    bvpdb = pg.DB('bvparam2012', 'localhost')

    # 1st coord sphere statistics
    files_1cs = {}
    print "from table...", table_1cs
    for an_name in an_list_1cs:
        an_list = []
        #output = stat_output(an_name, "all", "percentage >= 0.0", 100, \
        #                     table_1cs, stat_tab_1cs)
        #an_list.append((output[0], output[2]))
        output = stat_output(an_name, "top 50", "percentage >= " + threshold, \
                             100, table_1cs, stat_tab_1cs)
        an_list.append((output[0], output[2]))
        files_1cs[an_name] = an_list
   #out_order = ['atom_name', 'an_num', 'total', 'num_homo', 'median', 'mean', \
    out_order = ['atom_name', 'num_homo', 'total', 'median', 'mean', \
                 'stddev']
    stat_printer(files_1cs, out_order)

#    # w/in 4A statistics
#    files_4a = {}
#    print "from table...", table_4a
#    for an_name in an_list_4a:
#        an_list = []
#        #output = stat_output(an_name, "all", "percentage >= 0", 100, \
#        #                     table_4a, stat_tab_4a)
#        #an_list.append((output[0], output[2]))
#        output = stat_output(an_name, "top 50", "percentage >= "+threshold, \
#                             100, table_4a, stat_tab_4a)
#        an_list.append((output[0], output[2]))
#        files_4a[an_name] = an_list
#    stat_printer(files_4a, out_order)

#     # gnuplotting section
#    print "Generating gnuplot file..."
#    gnuplotter(files)
#    print "Running gnuplot..."
#    os.system('gnuplot ' + DIRECTORY + '{0}_histograms.gp'.format(atom_type))
    bvpdb.close()
    bvp.close()
