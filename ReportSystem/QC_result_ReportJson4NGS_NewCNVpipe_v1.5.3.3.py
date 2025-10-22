#!/usr/bin/env python
# -*- coding:utf-8 -*-

__date__ = u'2025/02/26 13:37'
__software__ = u'PyCharm'
__email__ = u'zhouyajun@yikongenomics.com'
__author__ = u'\u5468\u4e9a\u519b'
__describe__ = u'QC_result_ReportJson4NGS_NewCNVpipe_v1.5.3.3.py'

"""
# 集群脚本
# 仅适用于CNV1.6.6分析流程的对接，后续升级需要做局部调整
1）添加对NICS项目的支持，并且NICSSR项目不做AI评级，因此需要区分一下NICS项目的,具体控制见L378 -- 2021年12月15日
2）添加芯片和MT的质控模块 -- 2022年1月26日
3）针对CNV的结果进行优化，特别是结论解析需要读取_results.tsv -- 2022年2月17日
4）需要对每个样本添加CNV的结果注释信息 -- 2022年2月23日
    + 基本就是46,XN的临床建议: "该样本未见染色体拷贝数变异，请结合其他临床指征综合考虑。"
    + 非46,XN的临床建议: "该样本存在染色体拷贝数变异，可能导致发育异常，请结合其他临床指征综合考虑。"
    + NA的不写: -
    + 意外发现没有临床建议
5）针对MaReCs1和PGTSR项目的断点信息汇总输出，要求必须压缩包中包含断点信息的JSON，“_breakpoint_by_project.json” -- 2022年3月2日
6) 结果文件中添加一个对性别判断的CNV图谱，发生性染色体CNV异常则报告带性别的cnv全图，以及添加一个样本名称到结果报告中  -- 2022年3月2日
2022年3月15日 -- debug，发现存在："临床表现（HI/TS基因）": "经查DECIPHER、ClinGen、OMIM、DGV等数据库，未见剂量敏感基因相关疾病"的项目 ，L341
2022年3月31日 -- 下面是需要拓展的功能，具体要求如下：
7) result.json 添加CNV项目的各个样本的“染色体示意图”: ["XX_chr1_chr12_4M_mos.png", "XX_chr13_chr22_4M_mos.png"]
8) MaReCs2项目--即芯片MaReCs项目需要添加PGTAH的报告内容，单独优化PGTAH的脚本尝试嵌合到一个脚本中
9) 目前芯片POCX的项目已独立嵌入到分析流程中，可以直接读取压缩包中的XX_QC_Report.json，不需要单独再在调用一遍该脚本进行输出
10) 添加对QC数据包的识别功能，鉴于PGTM的数据逻辑进行包结果识别
2022年4月28日 -- 针对info.xlsx格式的特殊边换规则进行微调
芯片PGTM已升级添加了kinship、baf、pgtah，以及存在不分析的merge，当前脚本避免用于芯片的处理，仅处理NGS数据
2022年6月6日 -- 范闲CNV的结果存在换行符直接join字符串的情况
2022年6月7日 -- 为全部的注释添加句末标点符号，优化各个模块的功能类型，删除无意义的函数，函数归类处理，优化NGS逻辑模式以及添加送检JSON判断和启动流程控制逻辑
2022年6月10日 -- 调整NGS的PGTM（SNP）的QC输出逻辑
2022年9月22日 -- 适配1.6.9的CNV流程升级出现的报错， XX_MT.txt增加了一列信息"MT_CN_normalize"做个适配微调
2022年10月10日 -- 适配黎明的CNV分析流程，出具对应的QC和report的JSON文件，NICS需要保留功能使用之前的1.4.4的逻辑单独处理，其他胚胎类的全部产品都需要进行调整（这个调用使用同步调用脚本加以区分）
2022年12月13日 -- 针对/data/Project_2022Q4/Project_YKSZ_PGTM_221126_06F_01_0298,连锁输出三个文件夹的识别不的情况做出调整
2022年12月29日 -- /data/Project_2022Q4/Project_YKSZ_MaReCs1_221209_21F_01_0066 该项目的S_NAME转换纯数字类型前面含有0的情况下自动抹掉0，导致没有按字符串处理,"info.xlsx"文件按字符串输出即可
2023年1月30日 -- 关于PGTM的SNP逻辑升级，调整后所有的QC仅仅通过三个文件体现,之前的文件不在在汇总体现在excel文件中，
                具体如下: 1) *_QC_summary.report.txt-->给客户的质控文件;
                        2) *_QC_summary.txt-->给医学遗传部审核用的质控文件;
                        3) *_iSNP_summary.txt-->有效位点数统计文件
2023年7月13日 -- PGTM的conclusion_embryo文件输出到到QC_summary文件里面
2023年8月16日 -- CNV项目的LDPGTA结果需要补充到报告的JSON文件中
2023年11月6日 -- .chrom_llr.png --> .chrom_llr.upd.png 替换
2024年1月12日 -- _result_Report.json 每个样本增加 "variationResultInfoList": {"SangerAtlas": [], "variationDetectionResult":"" } SangerAtlas只用来增加CNV小图路径，另一个空值处理
2024年1月25日 -- 经过IT协调调整到_QC_Report.json，字段放入到ASA_QC_Data中，占用Key：VARIATION_RESULT_List: {{"pic_name": ""}}, 移除1月12日的添加内容
2024年11月27日 -- 针对_06AL_ 增加浅测序的_result_Report.json文件读取_total_cnv.json的对应分辨率的chimeDetectionResult这个值写入"CNV结果信息"
2025年2月26日 -- {"pipeline": "ld_pgta", "version": "v1.2.1.0"}升级，鉴于v1.2.0.0之后UPD同意更换为ROH
"""

import re
import sys
import json
import textwrap
import time
import socket
import argparse
import platform
import shutil
import logging
import zipfile
from pathlib import Path
from glob import glob
import subprocess
from operator import contains
from collections import defaultdict
import pandas as pd


def create_logging(log_name):
    logger_ = logging.getLogger(log_name)
    logger_.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s | L:%(lineno)d | %(message)s', '%Y-%m-%d %H:%M:%S')
    ch.setFormatter(formatter)
    logger_.addHandler(ch)
    return logger_


class NGS_SNP_QC(object):
    def __init__(self, zip_Path, PT):
        global logger
        self.zip_Path = zip_Path
        self.PT = PT
        self.Path_dir = Path(self.zip_Path).parent.absolute().as_posix().__str__()
        self.Pro_name = self.remove_suffix(Path(self.zip_Path).name, ".zip")
    
    @staticmethod
    def remove_suffix(input_string, suffix):
        # removesuffix只有python3.9才有，向下兼容构建当前函数
        if input_string.endswith(suffix):
            return input_string[:-len(suffix)]
        return input_string
    
    def GetSNP_QC_Files(self):
        dd_list = defaultdict(list)
        json_QC = {}
        with zipfile.ZipFile(self.zip_Path, 'r', compression=zipfile.ZIP_DEFLATED) as target:
            for each in target.namelist():
                # 注意NGS的SNP没有Kinship、BAF、PGTAH图
                if "_spgd_QC_data.json" in each:
                    with target.open(each) as f:
                        json_QC = json.loads(f.read())
                elif contains(each, "_vcfQC.txt") and contains(each, "data"):
                    dd_list["vcfQC"].append(each)
                elif contains(each, "_QC_summary.txt") and contains(each, "data"):
                    dd_list["QC_summary"].append(each)
                elif contains(each, "_QC_summary.report.txt") and contains(each, "data"):
                    dd_list["QC_summary.report"].append(each)
                elif contains(each, "_iSNP_summary.txt") and contains(each, "data"):
                    dd_list["iSNP_summary"].append(each)
                elif contains(each, ".jpg") and contains(each, "graph"):
                    dd_list["jpg"].append(each)
                elif contains(each, "_conclusion_embryo.txt"):
                    dd_list["conclusion_embryo"].append(each)
        new_dd = self.ConvertDD(dd_list)
        QC_Summary, xlsx_path = self.GetXlsxStatQC(new_dd)
        return json_QC, QC_Summary, xlsx_path, dd_list
    
    def ReadNGS_SNP_QC2NewNGS_QC_Report(self):
        # QC_input: Project_YKSZ_PGTM_210619_06EK_01_0644_spgd_QC_data.json
        QC_input, QC_Summary, xlsx_path, dd_list = self.GetSNP_QC_Files()
        dd_qc = {}
        res_SNP = {}
        zip_QC_summary_ = {}
        for x, y in QC_Summary.items():
            zip_QC_summary_[x] = Path(y).name
        if QC_input:
            # QC脚本出错就会缺图
            # SNP_graph = QC_input["report_graph"]
            SNP_graph = dd_list.get("jpg", [])
            for k, v in QC_input.items():
                if "Project_" in k:
                    dd_qc = QC_input[k]
            for k, v in dd_qc.items():
                # TODO 判断QC数据平台NGS
                if self.PT == "NGS" and v["platform"] != "NGS":
                    logger.warning("N)-p参数设置NGS与使用数据不对应，检查结果是否是NGS数据包！")
                    sys.exit(1)
                res_SNP[k] = {"NGS_QC_Data": self.SNP_data_format(v, k)}
            res_SNP["SNP连锁图谱"] = SNP_graph
            res_SNP["QC_summary"] = zip_QC_summary_
            return res_SNP, xlsx_path
        else:
            logger.warning("没有检查到_spgd_QC_data.json文件，可能是个CNV的项目")
            return {}, ""
    
    @staticmethod
    def SNP_data_format(demo, SN):
        """
        demo = {'platform': 'NGS', 'gene': 'HBA', 'mode': 'AR',
                'SNP_QC_result': {'raw_reads': '4,607,840',
                                  'GC%': '53',
                                  'high_quality_reads': '4,234,519',
                                  'high_quality_of_raw(%)': '91.90',
                                  'mapped_reads': '3,943,959',
                                  'mapping_rate(%)': '93.14',
                                  'mapped_of_raw(%)': '85.59',
                                  'detection_rate(%)': '57.961783439490446',
                                  'QC': 'PASS',
                                  'data_QC_information': '-',
                                  'on_target_rate': '-',
                                  'Ratio_Error': '65',
                                  'Count_Sites': '6',
                                  'Total_Error': '1',
                                  'Miss': '5',
                                  'ADO': '0',
                                  'Other': '0.9846',
                                  'Recombination(maternal)': '0.0769',
                                  'Recombination(paternal)': '0.0769'},
                'CNV_result': {'Karotype': '46,XN,-7p(p21.2→p11.2,~42Mb,×1,mos,~49%)',
                               'Sex_chromosome_karyotype': 'XX',
                               'Description': '嵌合疑似',
                               'QCresult': 'PASS',
                               'Sample_annotation_template1(cn)': '-'}, 'SNP_result': 'PASS'}
        """
        new = {"Sample_name": SN, "CNV图谱": []}
        for k, v in demo.items():
            if isinstance(v, dict):
                new.update(v)
            else:
                new[k] = v
        Gender = new.get("Sex_chromosome_karyotype", "")
        if Gender and "," not in Gender:
            new["Gender"] = Gender
        else:
            Gender1 = Gender.split(',')[0]
            new["Gender"] = Gender1
        new["high_quality_rate(%)"] = new.pop("high_quality_of_raw(%)")
        list_none = ["coverage_of_genome(%)", "valid_reads", "valid_reads_GC_content(%)",
                     "valid_reads_rate(%)", "CV(1000K_bin_size)"]
        for i in list_none:
            new[i] = ""
        new["data_QC_conclusion"] = {}
        new["data_QC_information"] = {}
        return new
    
    @staticmethod
    def ConvertDD(dd_list):
        main_name = set()
        temp = {}
        for each in dd_list["vcfQC"]:
            main_ = Path(each).name.replace("_vcfQC.txt", "")
            logger.info("女方_基因: {}".format(main_))
            main_name.add(main_)
        for i in main_name:
            temp[i] = defaultdict(list)
            for txt in dd_list.get("QC_summary", []):
                if contains(txt, i):
                    temp[i]["QC_summary"].append(txt)
            for txt2 in dd_list.get("QC_summary.report", []):
                if contains(txt2, i):
                    temp[i]["QC_summary.report"].append(txt2)
            for txt3 in dd_list.get("iSNP_summary", []):
                if contains(txt3, i):
                    temp[i]["iSNP_summary"].append(txt3)
            # 2023年7月13日 增加conclusion_embryo
            for txt_C in dd_list.get("conclusion_embryo", []):
                if contains(txt_C, i):
                    temp[i]["conclusion_embryo"].append(txt_C)
            for _jpg in dd_list.get("jpg", []):
                if contains(_jpg, i):
                    temp[i]["jpg"].append(_jpg)
        return temp
    
    def GetXlsxStatQC(self, new_dd):
        """
        {'Fengweimei_AR_HBA': defaultdict(<class 'list'>, {
        'QC_summary.report': ['data/Fengweimei_AR_HBA_QC_summary.report.txt'],
        'QC_summary': ['data/Fengweimei_AR_HBA_QC_summary.txt'],
        'iSNP_summary': ['data/Fengweimei_AR_HBA_iSNP_summary.txt'],
        'jpg': ['graph/Xiufeifei_XLR_AR.report.jpg', 'graph/Xiufeifei_XLR_AR.pedigree.jpg']})}
        输出QC_summary.xlsx文件
        """
        dd_QC = {}
        temp_list = []
        with zipfile.ZipFile(self.zip_Path, 'r', compression=zipfile.ZIP_DEFLATED) as target:
            for k, v in new_dd.items():
                xlsx_name = f"{self.Path_dir}/{self.Pro_name}_{k}_QC_summary.xlsx"
                xlWriter = pd.ExcelWriter(xlsx_name)
                try:
                    # sheet Title is more than 31 character
                    self.WriteModel(target, v, xlWriter, "QC_summary.report")  # 给客户的质控文件
                    self.WriteModel(target, v, xlWriter, "QC_summary")  # 给医学遗传部审核用的质控文件
                    self.WriteModel(target, v, xlWriter, "iSNP_summary")  # 有效位点数统计文件
                    self.WriteModel_conclusion_embryo(target, v, xlWriter, "conclusion_embryo")
                except Exception as E:
                    logger.error(E)
                # xlWriter.save()
                xlWriter.close()
                temp_list.append(xlsx_name)
                if Path(xlsx_name).exists():
                    dd_QC[k] = xlsx_name
        return dd_QC, temp_list
    
    def WriteModel(self, target, v, xlWriter, key_word="QC_summary"):
        # head_dd = {"QC_summary": 7,
        #            "QC_summary.report": 5,
        #            "iSNP_summary": 2}
        if v.get(f"{key_word}", []):
            logger.debug(f"检查连锁分析目录数量：{len(v.get(key_word, []))}")
            if len(v.get(f"{key_word}", [])) == 1:
                for i in v.get("{}".format(key_word), []):
                    head_num = self.TellHeadNumber(target, i)
                    logger.debug(f"{Path(i).as_posix()} 表头所在行数： {head_num}")
                    with target.open(i) as f:
                        # temp = pd.read_csv(f, header=head_dd[key_word], delimiter="\t")
                        temp = pd.read_csv(f, header=head_num - 1, delimiter="\t", dtype=str)
                        temp.to_excel(xlWriter, sheet_name="".join(f"{key_word}"[:31]), index=False)
            elif len(v.get(f"{key_word}", [])) >= 2:
                # TODO 遇到那种家系或者点位很多进行图片分割的情况使用sheet区分
                temp_df = []
                for i in v.get(f"{key_word}", []):
                    for number in range(1, len(v.get(f"{key_word}", [])) + 1, 1):
                        if re.findall(f"data{number}/|data_{number}/", i):
                            head_num = self.TellHeadNumber(target, i)
                            logger.debug(f"{Path(i).as_posix()} 表头所在行数： {head_num}")
                            with target.open(i) as f:
                                # 由于每个summary文件的表头都有说明，因此需要对表头文件单独处理，具体查看原文
                                # temp = pd.read_csv(f, header=head_dd[key_word], delimiter="\t")
                                temp = pd.read_csv(f, header=head_num - 1, delimiter="\t", dtype=str)
                                temp_df.append(temp)
                # 纵向合并每一行数据后删除非第一次的重复行，并输出一个新的df和新的index，排序按"#Name"列升序输出新df，保存excel
                new1 = pd.concat(temp_df).drop_duplicates(subset=None, keep="first", inplace=False, ignore_index=True)
                new = new1.sort_values(by=["#Name"], ascending=True, inplace=False)
                new.to_excel(xlWriter, sheet_name="".join(f"{key_word}"[:31]), index=False)
        else:
            logger.error("{}: 没有相关文件".format(key_word))
    
    @staticmethod
    def WriteModel_conclusion_embryo(target, v, xlWriter, key_word="conclusion_embryo"):
        if v.get(F"{key_word}", []):
            # 按女方和基因进行拆分了，文件通常只有一个，出现连锁拆分可能存在2个
            if len(v.get(f"{key_word}", [])) == 1:
                for i in v.get(f"{key_word}", []):
                    logger.debug(f"conclusion_embryo: {Path(i).as_posix()}")
                    file_name = Path(i).name.split(".txt")[0]
                    with target.open(i) as f:
                        temp = pd.read_csv(f, header=None, delimiter="\t", dtype=str, names=['#Name', 'conclusion_embryo'])
                        temp["file_name"] = file_name
                        temp.to_excel(xlWriter, sheet_name="".join(f"{key_word}"[:31]), index=False)
            elif len(v.get(F"{key_word}", [])) >= 2:
                # TODO 遇到那种家系或者点位很多按文件名字一致的基因进行合并,分别输出到单个基因的sheet中
                temp_df = []
                temp_dd = defaultdict(list)
                for i in v.get(f"{key_word}", []):
                    for number in range(1, len(v.get(f"{key_word}", [])) + 1, 1):
                        if re.findall(f"data{number}/|data_{number}/", i):
                            logger.debug(f"conclusion_embryo: {Path(i).as_posix()}")
                            file_name = Path(i).name.split(".txt")[0]
                            temp_dd[file_name].append(i)
                for k_FN, v_CE_list in temp_dd.items():
                    for each_CE in v_CE_list:
                        with target.open(each_CE) as f:
                            temp = pd.read_csv(f, header=None, delimiter="\t", dtype=str, names=['#Name', 'conclusion_embryo'])
                            temp["file_name"] = k_FN
                            temp_df.append(temp)
                # 纵向合并每一行数据后删除非第一次的重复行，并输出一个新的df和新的index，排序按"#Name"列升序输出新df，保存excel
                new1 = pd.concat(temp_df).drop_duplicates(subset=None, keep="first", inplace=False, ignore_index=True)
                new = new1.sort_values(by=["#Name"], ascending=True, inplace=False)
                new.to_excel(xlWriter, sheet_name="".join(f"{key_word}"[:31]), index=False)
        else:
            logger.error(f"{key_word}: 没有相关文件")
    
    @staticmethod
    def TellHeadNumber(target, file_byte):
        n = 0
        with target.open(file_byte) as f_byte:
            while True:
                line = f_byte.readline()
                n += 1
                if not line:
                    break
                if line.startswith(b"##"):
                    continue
                if line.startswith(b"#Name"):
                    return n


class CNV_QC(object):
    def __init__(self, QC_input, CNV_graph, json_CNV_info, MT_dd_list):
        self.QC_input = QC_input
        self.CNV_graph = CNV_graph
        self.json_CNV_info = json_CNV_info
        self.MT_dd_list = MT_dd_list
    
    @staticmethod
    def CNV_data_format(SN, Demo, CNV_graph, gender):
        """
        demo = {"Data_sts": {"CV(1000K_bin_size)": "0.138554728266295",
                             "coverage_of_genome(%)": "5.81973874857117",
                             "duplication_rate(%)": "3.7593",
                             "high_quality_rate(%)": "99.0670781890275",
                             "mapping_rate(%)": "98.7933027162738",
                             "raw_reads": "3015794",
                             "valid_reads": "2557001",
                             "valid_reads_GC_content(%)": "37.78743137142385",
                             "valid_reads_rate(%)": "84.7869914191752"},
                "data_QC_conclusion": {
                    "10M": "WARNING",
                    "1M": "FAIL",
                    "4M": "WARNING",
                    "Whole_Arm": "WARNING",
                    "Whole_Chromosome": "WARNING"},
                "data_QC_information": {
                    "10M":
                        "valid_reads_GC_content(%):INFO:Valid reads GC content is a little low(within normal range), be cautious for CNV of chr4,13,17,19,22,X;CV(1000K bin size):INFO:CV(1000K bin size) is a little high(within normal range), be cautious for CNV(segment/arm/chromsome) length less than or equal to 10Mb;",
                    "1M":
                        "valid_reads:WARNING:Insufficient valid reads, be cautious for CNV(segment/arm/chromsome) length less than or equal to 4Mb;valid_reads_GC_content(%):INFO:Valid reads GC content is a little low(within normal range), be cautious for CNV of chr4,13,17,19,22,X;CV(1000K bin size):INFO:CV(1000K bin size) is a little high(within normal range), be cautious for CNV(segment/arm/chromsome) length less than or equal to 10Mb;",
                    "4M":
                        "valid_reads_GC_content(%):INFO:Valid reads GC content is a little low(within normal range), be cautious for CNV of chr4,13,17,19,22,X;CV(1000K bin size):INFO:CV(1000K bin size) is a little high(within normal range), be cautious for CNV(segment/arm/chromsome) length less than or equal to 10Mb;",
                    "Whole_Arm":
                        "valid_reads_GC_content(%):INFO:Valid reads GC content is a little low(within normal range), be cautious for CNV of chr4,13,17,19,22,X;CV(1000K bin size):INFO:CV(1000K bin size) is a little high(within normal range), be cautious for CNV(segment/arm/chromsome) length less than or equal to 10Mb;",
                    "Whole_Chromosome":
                        "valid_reads_GC_content(%):INFO:Valid reads GC content is a little low(within normal range), be cautious for CNV of chr4,13,17,19,22,X;CV(1000K bin size):INFO:CV(1000K bin size) is a little high(within normal range), be cautious for CNV(segment/arm/chromsome) length less than or equal to 10Mb;"
                }}
        """
        # print(SN, gender)
        new = {"Sample_name": SN, "CNV图谱": [CNV_graph], "Gender": gender[SN]}
        for k, v in Demo.items():
            if k == "Data_sts":
                new.update(Demo[k])
            else:
                new[k] = v
        return new
    
    @staticmethod
    def Get_CNV_gender(json_CNV_info):
        dd = set()
        for k, v in json_CNV_info.items():
            gender = v.get("Sex_chromosome_karyotype", "")
            if gender == "" or gender == "N/A":
                continue
            else:
                dd.add(gender)
        if dd:
            if "," in list(dd)[0]:
                return list(dd)[0].split(",")[0]
            else:
                return list(dd)[0]
        else:
            return ""
    
    def ReadNGS_CVN_QC2NewNGS_QC_Report(self):
        # QC_input: Project_YKSZ_IBPGTM_211124_GW000023_01_0800/Project_YKSZ_IBPGTM_211124_GW000023_01_0800_total_data_sts.json
        # json_CNV_info : Project_YKSZ_IBPGTM_211124_GW000023_01_0800/Project_YKSZ_IBPGTM_211124_GW000023_01_0800_total_CNV.json
        # MT_dd_list: Project_YKSZ_PGTA_210513_07C_06_0601/local_report/Project_YKSZ_PGTA_210513_07C_06_0601_MT.txt
        gender_dd = {}
        for k_i, v_i in self.json_CNV_info.items():
            gender_dd[k_i] = self.Get_CNV_gender(v_i)
        # print(gender_dd)
        new = {}
        for k, v in self.QC_input.items():
            if self.MT_dd_list:
                # 可能不存在MT_QC，没有就空值
                new[k] = {"NGS_QC_Data": self.CNV_data_format(k, v, self.CNV_graph[k], gender_dd),
                          "MT_QC": self.MT_dd_list[k]}
            else:
                new[k] = {"NGS_QC_Data": self.CNV_data_format(k, v, self.CNV_graph[k], gender_dd),
                          "MT_QC": self.MT_dd_list}
        return new


class CNV_NewPipeReport(object):
    def __init__(self, ad: dict, TCD_json_CNV_info: dict, TCDIF_DF, NICS=None, bp=None, flag_CDR=False):
        self.ad = ad  # 补充信息 {"xx": {'样本名称': '53586_B1', 'CNV图谱': ['N/A'], '染色体示意图': ['N/A']}}
        self.NICS = NICS  # 默认不存在NICS项目，待后续升级逐步添加逻辑
        self.TCD_json_CNV_info = TCD_json_CNV_info  # _total_cnv.json 载入，根据分辨率和性别的验证选取对应的待报告的key值读取对应的注释信息
        self.TCDIF_DF = TCDIF_DF  # 意外发现的文件（“_Additional_cnv_details.tsv”）载入df
        self.bp = bp  # 针对PGTSR和MaReCs项目补充的断点信息追加，ONPGTSR也存在断点信息
        self.flag_CDR = flag_CDR  # chimeDetectionResult 获取06AL嵌合比例
        # self.CNV_Thumbnail = dict(CNV_Thumbnail)  # 为了生成这部分： "variationResultInfoList": {"SangerAtlas": [], "variationDetectionResult":"" }
    
    def ParseTCD_JSON2SampleDict(self):
        TCD = self.TCD_json_CNV_info.copy()
        info_ddd = {}
        for each_sample in TCD:
            sample_anno_info_temp = TCD[each_sample]
            res1 = sample_anno_info_temp.pop("Whole_Arm_mos_no_report_gender")
            res2 = sample_anno_info_temp.pop("Whole_Arm_mos_report_gender")
            sample_anno_info_temp.pop("Whole_Chromosome_mos_no_report_gender")
            sample_anno_info_temp.pop("Whole_Chromosome_mos_report_gender")
            # 过滤掉非嵌合 _nomos_ 保留含有嵌合的key
            sample_anno_info = {k: v for k, v in sample_anno_info_temp.items() if re.findall(r"_mos_", k)}
            x = ""
            if len(sample_anno_info) == 2:
                for i in sample_anno_info:
                    if '_no_report_gender' in i:
                        # 获取不带性别的指定分辨率下的注释信息
                        x = i
                logger.debug(x)
                res3 = sample_anno_info.pop(x)
                # 修改对应的key值,删除无效数据
                try:
                    res3.pop("Description_EN")
                    res3.pop("Ratings")
                    res3.pop("Morphology_Karyotype_Ratings")
                    res3.pop("Sample_annotation_template1_cn")
                    res3.pop("Sample_annotation_template1_en")
                    res3.pop("cnv_location")
                    res3.pop("report_graph")
                    res3.pop("chromosome_report_graph")
                    if self.flag_CDR:
                        res3["chimeDetectionResult"] = res3.get("chimeDetectionResult", "")
                    res3["CNV检测结果"] = res3.pop("Karotype")
                    res3["性别"] = res3.pop("Sex_chromosome_karyotype")
                    res3["结果解释"] = res3.pop("Description_CN")
                    res3["结果说明"] = self.Format_结果说明(res3.pop("result_description"))
                    res3["临床建议"] = self.Format_临床建议(res3.pop("Clinical_recommendations"))
                    res3["临床表现（HI/TS基因）"] = self.Format_剂量敏感基因(res3.pop("Dosage_Sensitivity_Gene_related_diseases"))
                    res3["临床表现（综合征）"] = self.Format_综合症(res3.pop("clinical_feature"))
                    res3["参考依据"] = self.FormatReference(res3.pop("Reference"))
                except Exception as ER:
                    logger.error(ER)
                res3.update(self.ad[each_sample])
                info_ddd[each_sample] = {"CNV结果信息": res3,
                                         "意外发现结果信息": self.ParseAdditional_cnv_details().get(each_sample, [])}
            elif len(sample_anno_info) == 0:
                # sample_anno_info["Whole_Arm_mos_no_report_gender"] = res1
                # sample_anno_info["Whole_Arm_mos_report_gender"] = res2
                # 修改对应的key值,删除无效数据
                try:
                    res1.pop("Description_EN")
                    res1.pop("Ratings")
                    res1.pop("Morphology_Karyotype_Ratings")
                    res1.pop("Sample_annotation_template1_cn")
                    res1.pop("Sample_annotation_template1_en")
                    res1.pop("cnv_location")
                    res1.pop("report_graph")
                    res1.pop("chromosome_report_graph")
                    if self.flag_CDR:
                        res1["chimeDetectionResult"] = res1.get("chimeDetectionResult", "")
                    res1["CNV检测结果"] = res1.pop("Karotype")
                    res1["性别"] = res1.pop("Sex_chromosome_karyotype")
                    res1["结果解释"] = res1.pop("Description_CN")
                    res1["结果说明"] = self.Format_结果说明(res1.pop("result_description"))
                    res1["临床建议"] = self.Format_临床建议(res1.pop("Clinical_recommendations"))
                    res1["临床表现（HI/TS基因）"] = self.Format_剂量敏感基因(res1.pop("Dosage_Sensitivity_Gene_related_diseases"))
                    res1["临床表现（综合征）"] = self.Format_综合症(res1.pop("clinical_feature"))
                    res1["参考依据"] = self.FormatReference(res1.pop("Reference"))
                except Exception as ER2:
                    logger.error(ER2)
                
                # info_ddd[each_sample] = res1
                res1.update(self.ad[each_sample])
                info_ddd[each_sample] = {"CNV结果信息": res1,
                                         "意外发现结果信息": self.ParseAdditional_cnv_details().get(each_sample, [])}
            elif len(sample_anno_info) > 2:
                # 兼容10M分辨率的项目
                x = ""
                for i in sample_anno_info.keys():
                    if "M_mos_no_report_gender" in i:
                        x = i
                logger.debug(x)
                res3 = sample_anno_info.pop(x)
                # 修改对应的key值,删除无效数据
                try:
                    res3.pop("Description_EN")
                    res3.pop("Ratings")
                    res3.pop("Morphology_Karyotype_Ratings")
                    res3.pop("Sample_annotation_template1_cn")
                    res3.pop("Sample_annotation_template1_en")
                    res3.pop("cnv_location")
                    res3.pop("report_graph")
                    res3.pop("chromosome_report_graph")
                    if self.flag_CDR:
                        res3["chimeDetectionResult"] = res3.get("chimeDetectionResult", "")
                    res3["CNV检测结果"] = res3.pop("Karotype")
                    res3["性别"] = res3.pop("Sex_chromosome_karyotype")
                    res3["结果解释"] = res3.pop("Description_CN")
                    res3["结果说明"] = self.Format_结果说明(res3.pop("result_description"))
                    res3["临床建议"] = self.Format_临床建议(res3.pop("Clinical_recommendations"))
                    res3["临床表现（HI/TS基因）"] = self.Format_剂量敏感基因(res3.pop("Dosage_Sensitivity_Gene_related_diseases"))
                    res3["临床表现（综合征）"] = self.Format_综合症(res3.pop("clinical_feature"))
                    res3["参考依据"] = self.FormatReference(res3.pop("Reference"))
                except Exception as ER:
                    logger.error(ER)
                res3.update(self.ad[each_sample])
                info_ddd[each_sample] = {"CNV结果信息": res3,
                                         "意外发现结果信息": self.ParseAdditional_cnv_details().get(each_sample, [])}
            else:
                logger.error(f"样本report字段长度异常：{len(each_sample)}")
                sys.exit(1)
        return info_ddd
    
    @staticmethod
    def FormatReference(Reference):
        Reference_str = ""
        if isinstance(Reference, list):
            Reference_str = "\n".join(Reference)
        elif isinstance(Reference, dict):
            for k, v in Reference.items():
                temp = [k]
                temp.extend(v)
                Reference_str += "{}\n\n".format("\n".join(temp))
        else:
            raise Exception("不识别的参考数据类型！")
        return Reference_str
    
    @staticmethod
    def Format_结果说明(result_description):
        RD = ""
        if isinstance(result_description, str):
            RD = result_description
        elif isinstance(result_description, list):
            seen = set()
            result = [] 
            for item in result_description:
                if item not in seen:
                    seen.add(item)
                    result.append(item)
            for i, v in enumerate(result, start=1):
                RD += f"{i}、{v}\n"
        else:
            raise Exception("结果说明既不是字符串格式也不是列表存储！")
        return RD
    
    @staticmethod
    def Format_临床建议(Clinical_recommendations):
        CR = ""
        if isinstance(Clinical_recommendations, str):
            CR = Clinical_recommendations
            # return CR
        else:
            raise Exception("临床建议的数据类型不识别！")
        return CR
    
    @staticmethod
    def Format_剂量敏感基因(Dosage_Sensitivity_Gene):
        # 2022年10月11日 -- 通过评审决议后续不再出具剂量敏感的相关注释
        # 结果永远是：["-"]
        return "\n".join(Dosage_Sensitivity_Gene)
    
    @staticmethod
    def Format_综合症(clinical_feature):
        CF = ""
        if isinstance(clinical_feature, list):
            if len(clinical_feature) > 1:
                for i, v in enumerate(clinical_feature, start=1):
                    CF += f"{i}、{v}\n"
            elif len(clinical_feature) == 1:
                if clinical_feature[0] in ["-", "/"]:
                    CF = "\n".join(clinical_feature)
            else:
                CF = "\n".join(clinical_feature)
        else:
            raise Exception("临床综合症的数据类型不识别！")
        return CF
    
    def Add_Breakpoint2SampleDict(self):
        sample_anno_dd = self.ParseTCD_JSON2SampleDict()
        # logger.debug(sample_anno_dd)
        if self.bp:
            sample_anno_dd["断点信息"] = self.bp
            return sample_anno_dd
        else:
            return sample_anno_dd
    
    def ParseAdditional_cnv_details(self):
        # 解析意外发现数据表到字典
        if self.TCDIF_DF.empty:
            return {}
        else:
            temp_dd = defaultdict(list)
            dict_TCDIF = self.TCDIF_DF.to_dict(orient="records")
            for i, each in enumerate(dict_TCDIF, start=1):
                temp_dd[each['Sample']].append({"序号": i,
                                                "提示CNV": each["Abnormal"],
                                                "染色体位置": "{0}:{1}-{2}".format(each['Chromosome'], each['Start_Pos'], each['End_Pos']),
                                                "临床建议": "",
                                                "临床表现（HI/TS基因）": "",
                                                "临床表现（综合征）": each["Syndrome(cn)"],
                                                "参考依据": "ClinGen：{0} \nDecipher：{1}".format(each["ClinGen"], each["Decipher"]),
                                                "结果说明": each["Phenotypes(cn)"]})
            return temp_dd


class Run_QC_Report(object):
    def __init__(self, input_zip, platform_type, analysis_type, zip_Path, project_name):
        self.input_zip = input_zip
        self.platform_type = platform_type
        self.analysis_type = analysis_type
        self.zip_Path = zip_Path
        self.project_name = project_name
    
    def Run(self):
        (input_zip, platform_type,
         analysis_type, zip_Path, project_name) = (self.input_zip, self.platform_type,
                                                   self.analysis_type, self.zip_Path, self.project_name)
        
        # PGD & IB_PGD
        PGD_project = ['19AFSLASR', '19ASLASR', '19AFMSLASR', 'PGTMFSPFREE', 'PGD', 'PGTM', 'PGTMSR', 'PGTMFFREE',
                       'PGTMFREE', 'PGDRD2000', 'PGTMRD2000', 'PGDM', 'PGTMM', 'PGTMMRD2000',
                       'NIPGTM', 'PGDC', 'PGTMC', 'PGD_JX', 'HLA-PGD', 'HLAPGTM', 'HLA-PGDF',
                       'HLAPGTMF', 'PGTMCRD2000', 'PGDF', 'PGDFRD2000', 'PGTMFRD2000', 'PGTMF', 'PGTMFT',
                       'PGDF-S', 'PGTMFS', 'PGDFSP', 'PGDSR', 'PGTMFSP', 'PGDA', 'PGTMA', 'PBDSR',
                       'IB-PGD', 'IB_PGD', 'IBPGTM', 'IB']
        n = 0
        snp, cnv, bp = "", "", ""
        TCD, TCDIF = "", ""  # TCD=total_CNV_detail， TCDIF=total_CNV_detail_IF
        SNP_QC_res, CNV_QC_res, MT_dd_list, json_BP_CNV = {}, {}, {}, {}
        json_QC_CNV, json_CNV_info = {}, {}
        withXY_mos_1M = {}
        noXY_mos_1M = {}
        CNV2LDPGTA = {}
        Chromosome_Diagram = defaultdict(list)
        UPD_Diagram = defaultdict(list)
        CNV_Thumbnail = defaultdict(list)
        if zipfile.is_zipfile(zip_Path):
            with zipfile.ZipFile(zip_Path, 'r', compression=zipfile.ZIP_DEFLATED) as target:
                for each in target.namelist():
                    if "_spgd_QC_data.json" in each:
                        snp = each
                        logger.info(F"检查到包含SNP的质控数据(*spgd_QC_data.json)! --> {each}")
                        # Demo : Project_YKSZ_PGTM_210619_06EK_01_0644_spgd_QC_data.json
                        n += 1
                    elif "{}_total_data_sts.json".format(project_name) in each:
                        cnv = each
                        logger.info(F"检查到包含CNV的质控数据(*total_data_sts.json)! --> {each}")
                        with target.open(each) as f_QC_CNV:
                            json_QC_CNV = json.loads(f_QC_CNV.read())
                        n += 1
                    elif "{}_total_cnv.json".format(project_name) in each:
                        logger.info(F"检查到包含total_cnv.json  --> {each}")
                        with target.open(each) as f_CNV:
                            json_CNV_info = json.loads(f_CNV.read())
                        n += 1
                    elif "_withXY_mos.png" in each:
                        # logger.info(F"检查到包含withXY_mos.png --> {each}")
                        key_name = Path(each).name.replace("_withXY_mos.png", "")
                        withXY_mos_1M[key_name] = each
                        n += 1
                    elif "_noXY_mos.png" in each:
                        # logger.info(F"检查到包含noXY_mos.png --> {each}")
                        key_name = Path(each).name.replace("_noXY_mos.png", "")
                        noXY_mos_1M[key_name] = each
                        n += 1
                    elif "_cytoband.png" in each:
                        # logger.info(F"检查到包含cytoband.png --> {each}")
                        key_name = Path(each).name.replace("_cytoband.png", "")
                        Chromosome_Diagram[key_name].append(each)
                    elif "info.xlsx" in each:
                        logger.info(F"检查到包含*info.xlsx --> {each}")
                        # TODO PGD&IB_PGD相关项目需要单独处理，其info表格只有3列，其他都是5列
                        #  -- 2022年4月28日 -- v9_bin/generate_info_and_fullgraph.py
                        pro_code = project_name.split("_")[1]
                        with target.open(each) as info_xlsx:
                            if pro_code in PGD_project:
                                try:
                                    info_sample = pd.read_excel(info_xlsx, header=None, na_values="-", dtype=str,
                                                                na_filter=False, keep_default_na=False,
                                                                names=["样本条码", "样本名称", "CNV结果"],
                                                                sheet_name='info', index_col="样本条码").to_dict(orient="index")
                                except Exception as ERROR:
                                    logger.error(F"*info.xlsx 异常强制使用5列模式读取文件！ -->{ERROR}")
                                    info_sample = pd.read_excel(info_xlsx, header=None, na_values="-", dtype=str,
                                                                na_filter=False, keep_default_na=False,
                                                                names=["样本条码", "样本名称", "性别", "CNV结果", "女方"],
                                                                sheet_name='info', index_col="样本条码").to_dict(orient="index")
                            else:
                                info_sample = pd.read_excel(info_xlsx, header=None, na_values="-", dtype=str,
                                                            na_filter=False, keep_default_na=False,
                                                            names=["样本条码", "样本名称", "性别", "CNV结果", "女方"],
                                                            sheet_name='info', index_col="样本条码").to_dict(orient="index")
                        n += 1
                    # TODO 注释读取_total_cnv.json和意外发现注释都直接读_Additional_cnv_details.tsv
                    elif "_Additional_cnv_details.tsv" in each:
                        logger.info(F"检查到包含*Additional_cnv_details.tsv --> {each}")
                        with target.open(each) as f2:
                            TCDIF = pd.read_csv(f2, header=0, delimiter="\t")
                        n += 1
                    elif "{}_data.tsv".format(project_name) in each:
                        logger.info(F"检查到包含*_data.tsv --> {each}")
                        n += 1
                        with target.open(each) as fMT:
                            MT_df = pd.read_csv(fMT, header=0, delimiter="\t", index_col="##SampleID", usecols=["##SampleID", "MT_CN"])
                            MT_df.columns = ["MT_copy_number"]
                            # ["Project_ID", "sample", "unique_mapped_reads_autosome", "unique_mapped_reads_MT", "MT_copy_number", "MT_CN_normalize"]
                        MT_dd_list = {}
                        res_MT = MT_df.to_dict(orient="index")
                        for i in res_MT:
                            each_sample_dd = res_MT[i]
                            each_sample_dd["Project_ID"] = project_name
                            each_sample_dd["sample"] = i
                            each_sample_dd["unique_mapped_reads_autosome"] = '-'
                            each_sample_dd["unique_mapped_reads_MT"] = '-'
                            each_sample_dd["MT_CN_normalize"] = '-'
                            MT_dd_list[i] = each_sample_dd
                    # TODO 增加断点信息整理输出的逻辑,最终断点信息输出在breakpoint_by_project.json -- “BreakPoint:{}”
                    elif "_BreakPoints_Kmeans.json" in each:
                        logger.info(F"检查到包含*_BreakPoints_Kmeans.json --> {each}")
                        if re.findall("_MaReCs1|_PGTSR", project_name):
                            with target.open(each) as BP_CNV:
                                json_BP_CNV = json.loads(BP_CNV.read())
                            pro, bp = self.parse_breakpoint_info(json_BP_CNV)
                        else:
                            bp = None
                        n += 1
                    elif "_cnv.LDPGTA.tsv" in each:
                        logger.info(F"检查到包含*_cnv.LDPGTA.tsv --> {each}")
                        n += 1
                        with target.open(each) as LDPGTA:
                            # LDPGTA_df = pd.read_csv(LDPGTA, delimiter="\t", header=0,
                            #                         dtype=str, na_values="", keep_default_na=False,
                            #                         usecols=["Sample", "LDPGTA_CNV_Type", "LDPGTA_UPD_arms"],
                            #                         index_col="Sample")
                            LDPGTA_df = pd.read_csv(LDPGTA, delimiter="\t", header=0,
                                                    dtype=str, na_values="", keep_default_na=False,
                                                    usecols=["Sample", "LDPGTA_CNV_Type", "LDPGTA_ROH_arms"],
                                                    index_col="Sample")
                        CNV2LDPGTA = LDPGTA_df.to_dict(orient="index")
                    elif ".chrom_llr.upd.png" in each:
                        n += 1
                        key_name = Path(each).name.replace(".chrom_llr.upd.png", "")
                        UPD_Diagram[key_name].append(each)
                    elif "other_check_graph" in each:
                        if not Path(each).name.endswith("_CNV.png"):
                            # if Path(each).name.endswith("_CNV.png"):
                            #     continue
                            barcode_name = Path(each).name.split("_chr")[0]  # _chr分割默认
                            CNV_Thumbnail[barcode_name].append(Path(each).as_posix())
        
        LDPGTA_dd_sample = self.Deal_LDPGTA(CNV2LDPGTA, dict(UPD_Diagram), dict(CNV_Thumbnail))
        logger.debug(LDPGTA_dd_sample)
        if cnv:
            additional_data = self.CNV_SupplementData(info_sample, withXY_mos_1M, noXY_mos_1M, Chromosome_Diagram)
        else:
            additional_data = {}
        if platform_type == "NGS":
            # TODO SNP的QC需要独立出一个类，QC_summary需要添加进去
            SNP_QC_res, QC_summary_path = NGS_SNP_QC(zip_Path, platform_type).ReadNGS_SNP_QC2NewNGS_QC_Report()
            # 已经添加XX_MT.txt,线粒体相关信息
            CNV_QC_res = CNV_QC(json_QC_CNV, withXY_mos_1M, json_CNV_info, MT_dd_list).ReadNGS_CVN_QC2NewNGS_QC_Report()
        else:
            # ASA 对PGTM|PGD仅输出_QC_Report.json, 对MaReCs|POCX|PGTAH单独写嵌入脚本
            # TODO 2022年5月9日 修订脚本，停止当前脚本对ASA芯片的支持
            logger.error("当前脚本只处理NGS的数据，ASA芯片的数据请使用ASA的脚本！")
            sys.exit(1)
        # logger.debug(CNV_QC_res.keys())
        if n > 0:
            # TODO Do Something
            TCD_json_CNV_info = json_CNV_info.copy()
            QC_file = self.DoSomething_PGTM(zip_Path, snp, cnv,
                                            SNP_QC_res, QC_summary_path, CNV_QC_res,
                                            platform_type, TCD_json_CNV_info, TCDIF,
                                            bp, additional_data, analysis_type,
                                            LDPGTA_dd_sample)
            with zipfile.ZipFile(input_zip, 'a', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as NT:
                if QC_file:
                    for i in QC_file:
                        for each in NT.namelist():
                            if Path(i).name in each:
                                # 已经存在的文件名无法覆盖，所以一旦存在已写入的文件需要删除后才能再次写入
                                logger.warning("存在待写入文件：%s -- 请尝试请使用zip -u更新压缩包！" % i)
                                sys.exit(1)
                                # cmd_zip = f"cd {Path(input_zip).parent.as_posix()} " \
                                #           f"&& zip -u {Path(input_zip).name} {Path(i).name}"
                                # ExecuteCMD(cmd_zip)
                        logger.info("开始写入数据：{}".format(i))
                        try:
                            NT.write(i, arcname=Path(i).name)
                        except Exception as e:
                            logger.error(e)
                else:
                    logger.warning("没找到待传递报告系统的JSON文件")
                    sys.exit(1)
        else:
            logger.warning("检查输入文件发现可能不是标准的zip结果压缩包，可能该项目仅含merge文件")
            sys.exit(1)
    
    @staticmethod
    def Deal_LDPGTA(CNV2LDPGTA: dict, UPD_Diagram: dict, CNV_Thumbnail: dict):
        # 增加ASA_QC_Data：{"亲本污染分析": {}}
        """
        {'01A230813PGTA01BTC01': {'LDPGTA_CNV_Type': 'Normal', 'LDPGTA_UPD_arms': '-'},
         '01A230813PGTA01BTC02': {'LDPGTA_CNV_Type': 'Normal', 'LDPGTA_UPD_arms': '-'},
         '01A230813PGTA01BTC03': {'LDPGTA_CNV_Type': 'Normal', 'LDPGTA_UPD_arms': '-'},
         '01A230813PGTA01BTC04': {'LDPGTA_CNV_Type': 'Normal', 'LDPGTA_UPD_arms': '-'},
         '01A230813PGTA01BTC05': {'LDPGTA_CNV_Type': 'Normal', 'LDPGTA_UPD_arms': '-'},
         '01A230813PGTA01BTC06': {'LDPGTA_CNV_Type': 'Normal', 'LDPGTA_UPD_arms': '-'},
         '01A230813PGTA01BTC07': {'LDPGTA_CNV_Type': 'Normal', 'LDPGTA_UPD_arms': '-'}}
         translation_dd = {"Normal": "二倍体样本，不是单倍体或者三倍体且没有明显的污染",
                          "UPD(Contam)": "有污染的UD",
                          "Contamination(mother)": "母源污染",
                          "Contamination(father)": "父源污染",
                          "Triploid(mother)": "母源三倍体",
                          "Triploid(father)":  "父源三倍体",
                          "Triploid": "三倍体",
                          "UPD(mother)": "母源UPD",
                          "UPD(father)": "父源UPD",
                          "UPD(mother,contam)": "有污染的母源UPD",
                          "UPD(father,contam)": "有污染的父源UPD"}
        Normal=未见异常
        Triploidy/Contamination=三倍体/污染（根据补测结果）
        Suspected UPD=全基因组UPD
        Suspected partial UPD=UPD(x)、UPD(xp/xq)  # 取LDPGTA_UPD_arms列结果 -- 亲本污染分析
        # 存在UPD需要输出XX.chrom_llr.png图片 -- BAF图谱
        """
        translate_LDPGTA_CNV_Type = {"Normal": "未见异常",
                                     "Triploidy/Contamination": "三倍体/污染（根据补测结果）",
                                     "Suspected UPD": "全基因组UPD",
                                     "Suspected partial UPD": "LDPGTA_UPD_arms",
                                     "partial UPD": "LDPGTA_UPD_arms",
                                     "Suspected ROH": "全基因组ROH",
                                     "Suspected partial ROH": "LDPGTA_UPD_arms",
                                     "partial ROH": "LDPGTA_UPD_arms",
                                     "Partial ROH": "LDPGTA_UPD_arms",
                                     "Contamination": "污染",
                                     "Triploidy": "三倍体",
                                     "": "-", "nan": "-", "UPD": "UPD"}
        # 2023年11月29日LD-pgta升级去掉Suspected字段
        res_dd = {}
        for k, v in CNV2LDPGTA.items():
            LDPGTA_CNV_Type = v.get("LDPGTA_CNV_Type", "")
            logger.debug((k, LDPGTA_CNV_Type))
            res = translate_LDPGTA_CNV_Type.get(str(LDPGTA_CNV_Type), "-")
            CNV_Thumbnail_pngs = CNV_Thumbnail.get(k, '')
            temp_CNV_info = defaultdict(list)
            if isinstance(CNV_Thumbnail_pngs, list):
                for each_png in CNV_Thumbnail_pngs:
                    temp_CNV_info["VARIATION_RESULT_List"].append({"pic_name": each_png,
                                                                   "Mutation_type": "",
                                                                   "geneName": "",
                                                                   "variationInfo": ""})
            elif isinstance(CNV_Thumbnail_pngs, str):
                temp_CNV_info["VARIATION_RESULT_List"].append({"pic_name": CNV_Thumbnail_pngs,
                                                               "Mutation_type": "",
                                                               "geneName": "",
                                                               "variationInfo": ""})
            else:
                logger.error(F"检查到不被识别的类型: {k} --> {CNV_Thumbnail_pngs}")
                sys.exit(1)
            
            if re.findall("LDPGTA_UPD_arms", res):
                # 此时选取LDPGTA_UPD_arms列的信息
                res2 = v.get("LDPGTA_ROH_arms", "-")
                res_dd[k] = {"ASA_QC_Data": {"BAF图谱": UPD_Diagram[k], "亲本污染分析": res2, "VARIATION_RESULT_List": dict(temp_CNV_info)["VARIATION_RESULT_List"]}}
            else:
                if re.findall("UPD", str(LDPGTA_CNV_Type)):
                    res_dd[k] = {"ASA_QC_Data": {"BAF图谱": UPD_Diagram[k], "亲本污染分析": res, "VARIATION_RESULT_List": dict(temp_CNV_info)["VARIATION_RESULT_List"]}}
                else:
                    # 无UPD不提供图片
                    res_dd[k] = {"ASA_QC_Data": {"BAF图谱": [], "亲本污染分析": res, "VARIATION_RESULT_List": dict(temp_CNV_info)["VARIATION_RESULT_List"]}}
        return res_dd
    
    @staticmethod
    def parse_breakpoint_info(json_BP_CNV: dict):
        global logger
        for pro, info in json_BP_CNV.items():
            logger.info("\n{}{}{}{}{}{}\n".format(pro, "--> 开始解析断点信息！", "\n罗氏易位：",
                                                  info["rob_common_breakpoint"], "\n断点信息：", info["samples_judged_breakpoint"]))
            if info["rob_common_breakpoint"] != "NA":
                # 罗氏易位的情况以罗氏断点为准
                bp = "\n".join(info["rob_common_breakpoint"].split("|"))
                return pro, bp
            else:
                # 非罗氏易位的情况结论直接取samples_judged_breakpoint
                bp = "\n".join(info["samples_judged_breakpoint"].split("|"))
                return pro, bp
    
    @staticmethod
    def CNV_SupplementData(info_sample: dict, withXY_mos_1M: dict, noXY_mos_1M: dict, Chromosome_Diagram: dict):
        """
        info_sample --> Project_XXinfo.xlsx
        {'24AR211103NICSA02B214': {'样本名称': '53586_B1', '性别': 'N/A', 'CNV结果': 'N/A', '女方': '李红婷'},
        '24AR211103NICSA02BC08': {'样本名称': '53586_1_1_1', '性别': 'XX,-Xq(×1,mos,~50%)', 'CNV结果': '46,XX,-Xq(×1,mos,~50%)', '女方': '李红婷'},
        '24AR211103NICSA02BC09': {'样本名称': '53586_1_2_1', '性别': 'XY,+X(×2,mos,~50%),-Y(×0,mos,~60%)', 'CNV结果': '46,XY,+X(×2,mos,~50%),-Y(×0,mos,~60%)', '女方': '李红婷'},
        '24AR211103NICSA02BC10': {'样本名称': '53586_1_3_1', '性别': 'XX', 'CNV结果': '46,XN', '女方': '李红婷'},
        '24AR211103NICSA02BC11': {'样本名称': '53586_1_4_1', '性别': 'XY', 'CNV结果': '46,XN', '女方': '李红婷'},
        '24AR211103NICSA02BC12': {'样本名称': '53586_1_5_1', '性别': 'XY', 'CNV结果': '46,XN,-14q(q21.1→q23.3,~28Mb,×1)', '女方': '李红婷'},
        '24AR211103NICSA02BC13': {'样本名称': '53586_1_7_1', '性别': 'XX', 'CNV结果': '46,XN', '女方': '李红婷'}}
        返回值仅保留待传递的数据
        """
        new_dd = {}
        for k, v in info_sample.items():
            CNV = v['CNV结果']
            # print(CNV)
            try:
                del v['性别']
                del v['女方']
            except Exception as Error_report:
                logger.error("{} -- 报告类型为PGD & IB_PGD的CNV项目文件中info.xlsx缺少一列信息，注意".format(Error_report))
                # continue
            logger.debug(f"{k}-->{CNV}")
            if CNV == "N/A":
                v["CNV图谱"] = noXY_mos_1M[k]
            else:
                # 目前的版本基本不用做判断是否存在XY异常，默认46.XN报告不体现性别
                # if CNV.find("Y", 6) != -1 or CNV.find("X", 6) != -1:
                if re.findall("XN", str(CNV)):
                    v["CNV图谱"] = noXY_mos_1M[k]
                elif re.findall("X|Y", str(CNV)):
                    # 从第七位检测到Y|X则报性别
                    v["CNV图谱"] = withXY_mos_1M[k]
                else:
                    v["CNV图谱"] = noXY_mos_1M[k]
            del v['CNV结果']
            # v["染色体示意图"] = Chromosome_Diagram[k]
            temp_list = sorted(Chromosome_Diagram[k], reverse=True)  # 逆序保证[1-12,13-22]
            # print(temp_list)
            v["染色体示意图"] = temp_list
            new_dd[k] = v
        # print(info_sample)
        return new_dd
    
    def DoSomething_PGTM(self, zip_file_path, SNP, CNV,
                         SNP_QC, QC_summary_path, CNV_QC_,
                         platform_type, TCD, TCDIF,
                         bp, ad, analysis_type,
                         LDPGTA_dd_sample):
        global logger
        project_name = Path(zip_file_path).name
        zip_file_path_name = Path(zip_file_path).as_posix().split('.zip')[0]
        if re.findall("_06AL_", project_name):
            logger.info(F"Code _06AL_ was recognized!")
            flag_chimeDetectionResult = True
        else:
            flag_chimeDetectionResult = False
        save_file1 = F"{zip_file_path_name}_QC_Report.json"
        save_file2 = F"{zip_file_path_name}_result_Report.json"
        if platform_type == "ASA":
            logger.error("仅含有ASA芯片质控数据,当前脚本不处理芯片数据！脚本退出！")
            sys.exit(1)
        elif platform_type == "NGS":
            # TODO 检查分析类型，CNV+SNP or SNP or CNV
            if SNP and not CNV:
                logger.info("仅含有SNP的质控数据")
                if analysis_type == "SNP":
                    with open(save_file1, 'w', encoding="utf-8") as F_SNP:
                        json.dump(SNP_QC, F_SNP, sort_keys=True, ensure_ascii=False, indent="\t")
                    res_list = [save_file1]
                    res_list.extend(QC_summary_path)
                    return res_list
                else:
                    logger.info("-a 参数值和包内容不符，脚本自动退出")
                    sys.exit(1)
            # 仅仅对含有CNV的项目提供输出小图的逻辑
            elif not SNP and CNV:
                logger.info("仅含有CNV的质控数据")
                if analysis_type == "CNV":
                    QC_dd = self.MergeDict(CNV_QC_, LDPGTA_dd_sample)
                    with open(save_file1, 'w', encoding="utf-8") as F_CNV:
                        json.dump(QC_dd, F_CNV, sort_keys=True, ensure_ascii=False, indent="\t")
                    
                    if re.findall("_MaReCs1|_PGTSR", project_name):
                        with open(save_file2, 'w', encoding="utf-8") as R_CNV:
                            json.dump(CNV_NewPipeReport(ad, TCD, TCDIF, NICS=None, bp=bp, flag_CDR=flag_chimeDetectionResult).Add_Breakpoint2SampleDict(),
                                      R_CNV, sort_keys=True, ensure_ascii=False, indent="\t")
                    else:
                        with open(save_file2, 'w', encoding="utf-8") as R_CNV:
                            json.dump(CNV_NewPipeReport(ad, TCD, TCDIF, NICS=None, flag_CDR=flag_chimeDetectionResult).Add_Breakpoint2SampleDict(),
                                      R_CNV, sort_keys=True, ensure_ascii=False, indent="\t")
                    return [save_file1, save_file2]
                else:
                    logger.error("-a 参数值和包内容不符，脚本自动退出")
                    sys.exit(1)
            elif SNP and CNV:
                logger.info("同时含有SNP和CNV的质控数据")
                if analysis_type in ['CNV_SNP', 'SNP_CNV']:
                    QC_dd = self.MergeDict(SNP_QC, CNV_QC_)
                    QC_dd2 = self.MergeDict(QC_dd, LDPGTA_dd_sample)
                    with open(save_file1, 'w', encoding="utf-8") as F_SC:
                        json.dump(QC_dd2, F_SC, sort_keys=True, ensure_ascii=False, indent="\t")
                    
                    if re.findall("_MaReCs1|_PGTSR", project_name):
                        with open(save_file2, 'w', encoding="utf-8") as R_CNV:
                            json.dump(CNV_NewPipeReport(ad, TCD, TCDIF, NICS=None, bp=bp, flag_CDR=flag_chimeDetectionResult).Add_Breakpoint2SampleDict(),
                                      R_CNV, sort_keys=True, ensure_ascii=False, indent="\t")
                    else:
                        with open(save_file2, 'w', encoding="utf-8") as R_CNV:
                            json.dump(CNV_NewPipeReport(ad, TCD, TCDIF, NICS=None, flag_CDR=flag_chimeDetectionResult).Add_Breakpoint2SampleDict(),
                                      R_CNV, sort_keys=True, ensure_ascii=False, indent="\t")
                    res_list = [save_file1, save_file2]
                    res_list.extend(QC_summary_path)
                    
                    return res_list
                else:
                    logger.error("-a 参数值和包内容不符，脚本自动退出")
                    sys.exit(1)
            else:
                logger.info("不包含SNP & CNV质控数据，稍后退出不作处理！")
                return []
    
    @staticmethod
    def MergeDict(a, b):
        import copy
        # 深度合并字典 a 和 b 的值
        merged_dict = copy.deepcopy(a)
        for key, value in b.items():
            if key in merged_dict and isinstance(merged_dict[key], dict) and isinstance(value, dict):
                merged_dict[key].update(value)
            else:
                merged_dict[key] = value
        return merged_dict


def arguments():
    parser = argparse.ArgumentParser(prog=__describe__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''
                                     \033[0;31;33m 脚本用于低深度NGS测序的PGD|PGT|MaReCs项目的QC&results文件整理合并  \033[0m'''))
    parser.add_argument('--version', action='version', version=__describe__)
    parser.add_argument('--ZipPath', '-z', dest='ZP', type=str, default=None, required=True,
                        help='ZIP file path')
    parser.add_argument('--PlatformType', '-p', dest='PT', type=str, default="NGS", choices=['NGS', 'ASA'],
                        help='NGS OR ASA')
    parser.add_argument('--AnalysisType', '-a', dest='AT', type=str, default="CNV",
                        choices=['SNP', 'CNV', 'CNV_SNP', 'SNP_CNV'], required=True,
                        help='Specify an analysis type CNV or SNP or CNV_SNP, default CNV Automatic report CNV project')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()


def main():
    global logger
    args = arguments()
    input_zip = args.ZP
    platform_type = args.PT
    analysis_type = args.AT  # SNP/CNV/SNP_CNV/CNV_SNP
    zip_Path = Path(input_zip).absolute()
    
    project_name = Path(zip_Path).name.split(".")[0].replace("Project_", "")
    # TODO 添加ban掉nics的逻辑，nics使用1.4.4版本出具报告
    if re.findall("ERT", project_name):
        logger.error("当前脚本不能处理ERT项目!")
        sys.exit(1)
    elif re.findall("POC|BLC|FCNV", project_name):
        logger.error("不能处理POC|BLC|FCNV！")
        sys.exit(1)
    elif re.findall("NICS|nics", project_name, flags=re.IGNORECASE):
        logger.error("新流程CNV不分析nics项目！")
        sys.exit(1)
    elif re.findall("RD|_RD_|_XSRD20_", project_name):
        logger.error("科研项目的CNV项目不处理！")
        sys.exit(1)
    elif re.findall("_Control_", project_name):
        logger.error("该流程Control项目直接不处理！")
        sys.exit(1)
    else:
        logger.info(f"识别到项目：{project_name}")
        # TODO 正式生产环境可以做微调，检查是否存在run.sh & analysis 目录判断分析版本1.6.X版本
        Run_QC_Report(input_zip, platform_type, analysis_type, zip_Path, project_name).Run()


if __name__ == "__main__":
    sys_platform = platform.system()
    HostName = socket.gethostname()
    format_time = time.strftime('%Y_%m_%d|%H:%M:%S', time.localtime(time.time()))
    logger = create_logging('LIMS2ReportSystemNGS')
    logger.info(f"Start: {__describe__}")
    main()
