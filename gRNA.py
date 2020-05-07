import re
import sys
from PyQt5 import QtWidgets, QtCore, QtGui
from fuzzysearch import find_near_matches


class UiGrna(object):
    def setupui(self, grna):
        grna.setObjectName("gRNA")
        grna.resize(650, 500)
        font = QtGui.QFont()
        font.setFamily("Adobe Arabic")
        font.setPointSize(7)
        grna.setFont(font)
        self.centralwidget = QtWidgets.QWidget(grna)
        self.centralwidget.setObjectName("centralwidget")
        self.gff_input = QtWidgets.QTextBrowser(self.centralwidget)
        self.gff_input.setGeometry(QtCore.QRect(280, 210, 251, 61))
        self.gff_input.setObjectName("gff_input")
        self.run_button = QtWidgets.QPushButton(self.centralwidget)
        self.run_button.setGeometry(QtCore.QRect(110, 390, 131, 61))
        font = QtGui.QFont()
        font.setFamily("Adobe Arabic")
        font.setPointSize(7)
        self.run_button.setFont(font)
        self.run_button.setObjectName("run_button")
        self.sel_gfffile_button = QtWidgets.QPushButton(self.centralwidget)
        self.sel_gfffile_button.setGeometry(QtCore.QRect(110, 210, 131, 61))
        font = QtGui.QFont()
        font.setFamily("Adobe Arabic")
        font.setPointSize(7)
        self.sel_gfffile_button.setFont(font)
        self.sel_gfffile_button.setObjectName("sel_gfffile_button")
        self.goalseq_input = QtWidgets.QTextBrowser(self.centralwidget)
        self.goalseq_input.setGeometry(QtCore.QRect(280, 40, 251, 61))
        self.goalseq_input.setObjectName("goalseq_input")
        self.run_out = QtWidgets.QTextBrowser(self.centralwidget)
        self.run_out.setGeometry(QtCore.QRect(280, 380, 251, 81))
        self.run_out.setObjectName("run_out")
        self.sel_goalseq_button = QtWidgets.QPushButton(self.centralwidget)
        self.sel_goalseq_button.setGeometry(QtCore.QRect(110, 40, 131, 61))
        font = QtGui.QFont()
        font.setFamily("Adobe 宋体 Std L")
        font.setPointSize(7)
        self.sel_goalseq_button.setFont(font)
        self.sel_goalseq_button.setObjectName("sel_goalseq_button")
        self.sel_simseq_button = QtWidgets.QPushButton(self.centralwidget)
        self.sel_simseq_button.setGeometry(QtCore.QRect(110, 125, 131, 61))
        font = QtGui.QFont()
        font.setFamily("Adobe Arabic")
        font.setPointSize(7)
        self.sel_simseq_button.setFont(font)
        self.sel_simseq_button.setObjectName("sel_simseq_button")
        self.simseq_input = QtWidgets.QTextBrowser(self.centralwidget)
        self.simseq_input.setGeometry(QtCore.QRect(280, 125, 251, 61))
        self.simseq_input.setObjectName("simseq_input")
        self.file_output_button = QtWidgets.QPushButton(self.centralwidget)
        self.file_output_button.setGeometry(QtCore.QRect(110, 295, 131, 61))
        font = QtGui.QFont()
        font.setFamily("Adobe Arabic")
        font.setPointSize(7)
        self.file_output_button.setFont(font)
        self.file_output_button.setObjectName("file_output_button")
        self.outfile_choose = QtWidgets.QTextBrowser(self.centralwidget)
        self.outfile_choose.setGeometry(QtCore.QRect(280, 295, 251, 61))
        self.outfile_choose.setObjectName("outfile_choose")
        grna.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(grna)
        self.menubar.setGeometry(QtCore.QRect(50, 30, 577, 25))
        self.menubar.setObjectName("menubar")
        grna.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(grna)
        self.statusbar.setObjectName("statusbar")
        grna.setStatusBar(self.statusbar)

        self.retranslateui(grna)
        QtCore.QMetaObject.connectSlotsByName(grna)

    def retranslateui(self, grna):
        _translate = QtCore.QCoreApplication.translate
        grna.setWindowTitle(_translate("gRNA", "基因编辑gRNA设计软件"))
        self.run_button.setText(_translate("gRNA", "开始设计"))
        self.sel_gfffile_button.setText(_translate("gRNA", "目标序列注释"))
        self.sel_goalseq_button.setText(_translate("gRNA", "目标序列文件"))
        self.sel_simseq_button.setText(_translate("gRNA", "相近序列文件"))
        self.file_output_button.setText(_translate("gRNA", "选择输出文件"))


class MainUi(QtWidgets.QMainWindow, UiGrna):
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        UiGrna.__init__(self)
        self.setupui(self)
        self.goalfile_input = ""
        self.fasta_input = ""
        self.gfffile_input = ""
        self.file_output = ""
        self.sel_goalseq_button.clicked.connect(self.selgoal)
        self.sel_simseq_button.clicked.connect(self.selsimseq)
        self.sel_gfffile_button.clicked.connect(self.selanno)
        self.file_output_button.clicked.connect(self.selout)
        self.run_button.clicked.connect(self.run)
        self.fileDialog = QtWidgets.QFileDialog(self)

    def selgoal(self):
        self.goalfile_input, _filter = self.fileDialog.getOpenFileName()
        self.goalseq_input.setText("".join(self.goalfile_input) + "\n")
        print("input"+"".join(self.goalfile_input)+"\n")

    def selsimseq(self):
        self.fasta_input, _filter = self.fileDialog.getOpenFileName()
        self.simseq_input.setText("".join(self.fasta_input) + "\n")
        print("input"+self.fasta_input+"\n")

    def selanno(self):
        self.gfffile_input, _filter = self.fileDialog.getOpenFileName()
        self.gff_input.setText("".join(self.gfffile_input) + "\n")
        print("input"+self.gfffile_input+"\n")

    def selout(self):
        self.file_output, _filter = self.fileDialog.getSaveFileName()
        self.outfile_choose.setText("".join(self.file_output)+"\n")
        print("output"+self.file_output+"\n")

    def run(self):
        out = "Seq Searching Job Starts"+"\n"
        self.run_out.setText(out)
        self.run_out.repaint()
        AR = AnnoResult
        grna_list = AR.read_file(self.goalfile_input)
        seq_gene = AR.readFasta(self.fasta_input)
        seq_list = []
        for i in seq_gene:
            seq_list.append(i)
        out += "Seq Scoring Job Starts"+"\n"
        self.run_out.setText(out)
        self.run_out.repaint()
        match_list = AR.matchesRNA(seq_list, grna_list, self.goalfile_input)
        site_list = AR.seq2mismatch(match_list)
        site_score = AR.sites2score(site_list)
        score_grade_output = AR.end_score(site_score)
        score_list = score_grade_output[0]
        grade_list = score_grade_output[1]
        dict_score = dict(score_list)
        dict_grade = dict(grade_list)
        out += "Seq Annotating Job Starts"+"\n"
        self.run_out.setText(out)
        self.run_out.repaint()
        dict_gff = AR.readgff(self.gfffile_input)
        anno_output = AR.exon_pos(dict_gff, score_grade_output)
        if anno_output==None:
            dict_anno = None
        else:
            dict_anno = dict(anno_output)
        AR.write_file(self.file_output, grna_list, dict_score, dict_grade, dict_anno)
        out += "All Job Done"+"\n"
        self.run_out.setText(out)
        self.run_out.repaint()
        print("All Job Done.")


class AnnoResult:
    @staticmethod
    def read_file(file_in):
        with open(file_in, "r") as f:
            seq = ""
            lines = f.readlines()
            for line in lines:
                if not line.startswith(">"):
                    seq += line.strip()
            seq = seq.upper()
            gg_site_list = []
            cc_site_list = []
            gg_starts = []
            cc_starts = []
            for g in re.compile("(?=GG)").finditer(seq):
                gg_start = g.start()
                gg_starts.append(gg_start)
            for gg_start in gg_starts:
                if gg_start >= 21:
                    output_line = seq[(gg_start - 21):(gg_start + 2)]
                    if "N" not in output_line:
                        gg_site_list.append([gg_start - 21, output_line])
            for c in re.compile("(?=CC)").finditer(seq):
                cc_start = c.start()
                cc_starts.append(cc_start)
            for cc_start in cc_starts:
                if len(seq) - cc_start < 21:
                    pass
                else:
                    output_line = seq[(cc_start):(cc_start + 23)]
                    if "N"not in output_line:
                        cc_site_list.append([cc_start, output_line])
        return gg_site_list, cc_site_list

    @staticmethod
    def readFasta(f):
        with open(f, 'r') as FA:
            seqName, seq = '', ''
            while 1:
                line = FA.readline()
                line = line.strip('\n')
                if (line.startswith('>') or not line) and seqName:
                    yield ((seqName, seq))
                if line.startswith('>'):
                    seqName = line[1:]
                    seq = ''
                else:
                    seq += line
                if not line:
                    break

    @staticmethod
    def matchesRNA(sq, l, file):
        def complement(letters1):
            letters2=[]
            basecomplemt = {"A": "T", "T": "A", "G": "C", "C": "G", "a": "t", "t": "a", "g": "c", "c": "g"}
            for base in letters1:
                if base in basecomplemt:
                    base = basecomplemt[base]
                else:
                    base = base
                letters2.append(base)
            return "".join(letters2)

        def revcomp(se):
            return complement(se)[::-1]

        with open(file, "r") as f:
            seq = ""
            lines = f.readlines()
            for line in lines:
                if not line.startswith(">"):
                    seq += line.strip()
            targ_seq = seq.upper()
            ant_targ_seq = revcomp(targ_seq)
            seq_len = len(targ_seq)
        gg_fuzzysearch = []
        cc_fuzzysearch = []
        gg_list = l[0]
        cc_list = l[1]
        for seq in gg_list:
            patten_seq = seq[1]
            find_string = []
            ant_patten_seq = revcomp(patten_seq)
            for item in sq:
                dna_seq = item[1]
                try:
                    start = dna_seq.index(targ_seq or ant_targ_seq)
                    start_int = int(start)
                    end = start + seq_len - 1
                    pre_seq = dna_seq[0:start_int]
                    mid_seq0 = dna_seq[start_int:(end + 1)]
                    mid_seq = "N" * seq_len
                    pro_seq = dna_seq[(end + 1):]
                    dna_seq = pre_seq + mid_seq + pro_seq
                except:
                    dna_seq = item[1]
                matches = find_near_matches(patten_seq, dna_seq, max_substitutions=4, max_insertions=0, max_deletions=0)
                for c in matches:
                    find_string.append(dna_seq[c.start:c.end])
                ant_matches = find_near_matches(ant_patten_seq, dna_seq, max_substitutions=4, max_insertions=0,
                                                max_deletions=0)
                for c in ant_matches:
                    ant_find_seq = dna_seq[c.start:c.end]
                    find_seq = revcomp(ant_find_seq)
                    find_string.append(find_seq)
            if find_string == []:
                find_string = "100.0"
            gg_fuzzysearch.append([seq, find_string])
        for seq in cc_list:
            find_string = []
            patten_seq = seq[1]
            ant_patten_seq = revcomp(patten_seq)
            for item in sq:
                dna_seq = item[1]
                matches = find_near_matches(ant_patten_seq, dna_seq, max_substitutions=4, max_insertions=0,
                                            max_deletions=0)
                for c in matches:
                    find_string.append(dna_seq[c.start:c.end])
                ant_matches = find_near_matches(patten_seq, dna_seq, max_substitutions=4, max_insertions=0,
                                                max_deletions=0)
                for c in ant_matches:
                    ant_find_seq = dna_seq[c.start:c.end ]
                    find_seq = revcomp(ant_find_seq)
                    find_string.append(find_seq)
            if find_string==[]:
                find_string = "100.0"
            cc_fuzzysearch.append([seq, [ant_patten_seq, find_string]])
        return gg_fuzzysearch, cc_fuzzysearch

    @staticmethod
    def seq2mismatch(L):
        seq_out = []
        gg_list = L[0]
        cc_list = L[1]
        for seq in gg_list:
            site_out = []
            site_query_dna = seq[0]
            query_dna = seq[0][1]
            ref_dnas = seq[1]
            if ref_dnas == "100.0":
                seq_out.append(seq)
            else:
                for sim_seq in ref_dnas:
                    site = []
                    ref_dna = sim_seq
                    if len(query_dna) == len(ref_dna):
                        i = 0
                        while i < len(query_dna):
                            if query_dna[i] != ref_dna[i]:
                                if i < 20:
                                    site.append(i)
                            i = i + 1
                    site_out.append(site)
                seq_out.append([site_query_dna, site_out])
        for seq in cc_list:
            site_out = []
            site_query_dna = seq[0]
            site = site_query_dna[0]
            query_dna = seq[1][0]
            site_ant_dna = tuple([site, query_dna])
            ref_dnas = seq[1][1]
            if ref_dnas == "100.0":
                seq_out.append([site_ant_dna, ref_dnas])
            else:
                for sim_seq in ref_dnas:
                    site = []
                    ref_dna = sim_seq
                    if len(query_dna) == len(ref_dna):
                        i = 0
                        while i < len(query_dna):
                            if query_dna[i] != ref_dna[i]:
                                if i < 20:
                                    site.append(i)
                            i = i + 1
                    site_out.append(site)
                seq_out.append([site_ant_dna, site_out])
        return seq_out

    @staticmethod
    def sites2score(mis_sites):
        seq_score = []
        for item in mis_sites:
            site_name = item[0]
            all_mis_site = item[1]
            if all_mis_site == "100.0":
                seq_score.append(item)
            else:
                score = []
                for site in all_mis_site:
                    mis_site = site
                    weight_all = (0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732,
                                  0.828, 0.615, 0.804,0.685, 0.583)
                    mis_site_num = len(mis_site)
                    score_p3 = 1
                    if mis_site_num >= 1:
                        score_p3 = 1 / (mis_site_num * mis_site_num)
                    score_p2 = 1
                    if mis_site_num >= 2:
                        distance = []
                        i = 0
                        while i < mis_site_num - 1:
                            dis = mis_site[i + 1] - mis_site[i]
                            distance.append(dis)
                            i = i + 1
                        sum_dis = 0
                        a = 0
                        while a < len(distance):
                            sum_dis = sum_dis + distance[a]
                            a = a + 1
                        mean_dis = sum_dis / len(distance)
                        score_p2 = 1 / ((((19 - mean_dis) / 19) * 4) + 1)
                    score_p1 = 1
                    b = 0
                    while b < mis_site_num:
                        score_p1 = score_p1 * (1 - weight_all[mis_site[b]])
                        b = b + 1
                    score_hit = score_p1 * score_p2 * score_p3 * 100
                    score.append(score_hit)
                seq_score.append((site_name, score))
        return seq_score

    @staticmethod
    def end_score(dif_list):
        sum_score = 0
        grade_list = []
        score_list = []
        for c in dif_list:
            gRNA = tuple(c[0])
            scores = c[1]
            if scores == "100.0":
                grade = "green"
                score_list.append((gRNA, scores))
                grade_list.append((gRNA, grade))
            else:
                for item in scores:
                    sum_score += item
                grna_score = 100/(100+sum_score)
                score_list.append((gRNA, grna_score))
                if grna_score < 0.5:
                    grade = "red"
                else:
                    grade = "green"
                grade_list.append((gRNA, grade))
        score_list = sorted(score_list, key=lambda d: d[0][0])
        grade_list = sorted(grade_list, key=lambda d: d[0][0])
        return score_list, grade_list

    @staticmethod
    def readgff(fg):
        with open(fg, "r")as gffin:
            gene1gff = {}
            if not fg == None:
                fg = gffin.readlines()
                d_keys = []
                d_mrna = []
                exon_starts = []
                exons = []
                exons1 = []
                for line in fg:
                    if not line.startswith("#"):
                        lin = line.rstrip().split("\t")
                        type1 = lin[2]
                        start = lin[3]
                        attribute = lin[8]
                        if type1 == "gene":
                            hit1 = re.search(r"ID\=([^;]+)", attribute)
                            dna = hit1.group(1)
                            d_keys.append(dna)
                        elif type1 == "mRNA":
                            hit2 = re.search(r"ID\=([^;]+)",attribute)
                            mrna = hit2.group(1)
                            d_mrna.append(mrna)
                        elif type1 == "exon":
                            exons.append(lin)
                        elif type1 == "three_prime_UTR" or "five_prime_UTR" or "CDS":
                            exon_starts.append(start)
                            exons.append(lin)
                        else:
                            pass
                for line in exons:
                    start = line[3]
                    type = line[2]
                    try:
                        if type == "exon" and start in exon_starts:
                            pass
                        else:
                            exons1.append(line)
                    except:
                        exons1.append(line)
                for i in range(len(d_mrna)):
                    valuei = []
                    for item in exons1:
                        if d_mrna[i] in item[8]:
                            valuei.append(item)
                    gene1gff[d_keys[i]] = valuei
            return gene1gff

    @staticmethod
    def exon_pos(d, li):
        output_list = []
        if d==None :
            output_list = None
        else:
            seq_len = int(23)
            site_list1 = []
            gg_list = li[0]
            cc_list = li[1]
            for item in gg_list:
                site_list1.append(item[0])
            for item in cc_list:
                site_list1.append(item[0])
            site_list1 = sorted(site_list1, key=lambda a: a[0])
            for key in d.keys():
                site_names = []
                exon_starts = []
                exon_ends = []
                utr_starts = []
                utr_ends =[]
                for value in d[key]:
                    type = value[2]
                    if type == "three_prime_UTR":
                        utr_starts.append(int(value[3]))
                        utr_ends.append(int(value[4]))
                    elif type == "five_prime_UTR":
                        utr_starts.append(int(value[3]))
                        utr_ends.append(int(value[4]))
                    else:
                        exon_starts.append(int(value[3]))
                        exon_ends.append(int(value[4]))
                exon_starts = sorted(exon_starts, reverse=True)
                exon_ends = sorted(exon_ends, reverse=True)
                for i in range(len(site_list1)):
                    site_name = tuple(site_list1[i])
                    site_names.append(site_name)
                    item = site_list1[i][0]
                    for c in range(len(exon_starts)):
                        try:
                            if item + seq_len <= utr_ends[0]:
                                output_list.append((site_name, "in the UTR of gene%s" % key))
                                continue
                            elif item>=utr_starts[1]:
                                output_list.append((site_name, "in the UTR of gene%s" % key))
                                continue
                            elif item <= utr_starts[1] and item+seq_len >=utr_starts[1]:
                                output_list.append((site_name, "between CDS and UTR of %s" % key))
                                continue
                            elif exon_starts[c] <= item <= exon_ends[c]:
                                if item + seq_len >= exon_ends[c]:
                                    output_list.append((site_name, "between CDS %d and intron of gene%s" % ((len(exon_starts) -
                                                                                                             c), key)))
                                    continue
                                else:
                                    output_list.append((site_name, "in the CDS %d of gene%s" % ((len(exon_starts) - c), key)))
                                    continue
                            elif item <= exon_starts[c] and item + seq_len >= exon_starts[c]:
                                if item <= utr_ends[0] and item + seq_len >= utr_ends[0]:
                                    output_list.append((site_name, "between UTR and CDS of gene%s" % key))
                                    continue
                                else:
                                    output_list.append((site_name, "between intron and CDS %d of gene%s" % ((len(exon_starts) - c),
                                                                                                         key)))
                                    continue
                            else:
                                continue
                        except:
                            if exon_starts[c] <= item <= exon_ends[c]:
                                if item + seq_len >= exon_ends[c]:
                                    output_list.append((site_name, "between exon %d and intron of gene%s" % ((len(exon_starts) -
                                                                                                          c), key)))
                                    continue
                                else:
                                    output_list.append((site_name, "in the exon %d of gene%s" % ((len(exon_starts) - c), key)))
                                    continue
                            elif item <= exon_starts[c] and item + seq_len >= exon_starts[c]:
                                output_list.append((site_name, "between intron and exon %d of gene%s" % ((len(exon_starts) - c),
                                                                                                 key)))
                                continue
                            else:
                                continue
                anno_site=[]
                for line in output_list:
                    anno_site.append(line[0])
                for name in site_names:
                    if not name in anno_site:
                        output_list.append((name, "in the intron of gene%s" % key))
        return output_list

    @staticmethod
    def write_file(out, lis, dict1, dict2, dict3):
        with open(out, "w")as output_file:
            output_file.write("number\tseq\tstart position\tstrand\tannotation\tscore\tgrade\n")
            list=[]
            Keys = []
            gg_site = []
            for item in lis[0]:
                gg_site.append(item[0])
            for i in dict1.keys():
                Keys.append(i)
            for i in range(len(dict1)):
                site_list = Keys[i]
                seq = site_list[1]
                site = site_list[0]
                if site in gg_site:
                    strand = "+"
                else:
                    strand = "-"
                if dict3==None:
                    ann = "NA"
                else:
                    if site_list in dict3.keys():
                        ann = dict3[site_list]
                    else:
                        ann = "."
                score = dict1[site_list]
                grade = dict2[site_list]
                c = (str(i+1), seq, str(site), strand, ann, str(score), grade)
                list.append(c)
            list1 = sorted(list,key=lambda e:e[5],reverse=True)
            for item in list1:
                item = "\t".join(item)
                output_file.write(item+"\n")
            output_file.close()


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = MainUi()
    window.show()
    sys.exit(app.exec_())
