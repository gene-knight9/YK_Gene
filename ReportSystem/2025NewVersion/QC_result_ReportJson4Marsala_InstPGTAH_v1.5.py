#!/usr/bin/env python
# -*- coding:utf-8 -*-

__date__ = u'2024/12/10 13:37'
__software__ = u'PyCharm'
__email__ = u'zhouyajun@yikongenomics.com'
__author__ = u'\u5468\u4e9a\u519b'
__describe__ = u'QC_result_ReportJson4Marsala_InstPGTAH_v1.3.py'

# TODO
#  1) 继承QC_result_ReportJson4NGS_NewCNVpipe_v1.5.1.py脚本的绝大多数CNV模块的功能;
#  2) 由于不涉及8基因SNP的分析模块,该板块冗余直接被删除;
#  3) 目前仅仅支持Inst-PGTA的Marsala项目分析结果,即无连锁分析的项目使用
#  4) 2023年10月13日 -- 由于PGTA合并到marsala流程，导致部分功能变动较大，进行模块功能微调
#  5) 2023年12月8日 -- 根据报告端反馈适当调整QC_report.JSON对接输入文件,兼容10月后的输出文件版本
#  6) 2024年2月1日 -- 经过IT协调调整到_QC_Report.json，字段放入到ASA_QC_Data中，占用Key：VARIATION_RESULT_List: {{"pic_name": ""}}
# 7) 2025-06-19 --增加"CNV检测结果（来源）"到_result_Report.json中

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
import colorlog
import pandas as pd
import os

def create_logging(log_name):
    logger_ = logging.getLogger(log_name)
    logger_.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = colorlog.ColoredFormatter(
        '%(log_color)s%(asctime)s | %(name)s | L:%(lineno)d | %(levelname)-8s%(reset)s %(message)s', '%Y-%m-%d %H:%M:%S',
        log_colors={
            'DEBUG': 'cyan',
            'INFO': 'green',
            'WARNING': 'yellow',
            'ERROR': 'red',
            'CRITICAL': 'black,bg_white'},
        reset=True,
        style='%')
    ch.setFormatter(formatter)
    logger_.addHandler(ch)
    return logger_


def ExecuteCMD(CMD):
    try:
        res = subprocess.check_call(CMD, shell=True)
        return res
    except subprocess.CalledProcessError as exc:
        logger.error('ReturnCode: {}'.format(exc.returncode))
        logger.error('CMD: {}'.format(exc.cmd))
        logger.error('Output: {}'.format(exc.output))
        return exc.returncode


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
    def __init__(self, ad: dict, TCD_json_CNV_info: dict, TCDIF_DF, temp_roh, NICS=None, bp=None, flag_CDR=False,orgin_dict = {}):
        self.ad = ad  # 补充信息 {"xx": {'样本名称': '53586_B1', 'CNV图谱': ['N/A'], '染色体示意图': ['N/A']}}
        # self.res_DF = res_DF  # {'07C210513PGTA06CE02': {'CNV检测结果': '46,XN', '结果解释': '未见染色体拷贝数异常'}} 数据来自：_cnv.tsv
        self.NICS = NICS  # 默认不存在NICS项目，待后续升级逐步添加逻辑
        self.TCD_json_CNV_info = TCD_json_CNV_info  # _total_cnv.json 载入，根据分辨率和性别的验证选取对应的待报告的key值读取回去对应的注释信息
        self.TCDIF_DF = TCDIF_DF  # 意外发现的文件（“_Additional_cnv_details.tsv”）载入df
        self.bp = bp  # 针对PGTSR和MaReCs项目补充的断点信息追加，ONPGTSR也存在断点信息
        self.flag_CDR = flag_CDR  # chimeDetectionResult 获取06AL嵌合比例
        self.temp_roh = temp_roh
        self.origin_dict = orgin_dict
    

    
    @staticmethod
    def add_snp_roh(low_res, roh_str):
        new_roh = roh_str.split("\n")
        if isinstance(low_res, str):
            return new_roh
        elif isinstance(low_res, list):
            low_res.extend(new_roh)
            return low_res
        else:
            return low_res
    
    
    def ParseTCD_JSON2SampleDict(self):
        TCD = self.TCD_json_CNV_info.copy()
        info_ddd = {}
        for each_sample in TCD:
            sample_anno_info_temp = TCD[each_sample]
            res1 = sample_anno_info_temp.pop("Whole_Arm_mos_no_report_gender")
            res2 = sample_anno_info_temp.pop("Whole_Arm_mos_report_gender")
            sample_anno_info_temp.pop("Whole_Chromosome_mos_no_report_gender")
            sample_anno_info_temp.pop("Whole_Chromosome_mos_report_gender")
            sample_orgin_result = None
            if each_sample in self.origin_dict.keys():
                sample_orgin_result = []
                sample_orgin_dict = self.origin_dict[each_sample]
                for chrom in sample_orgin_dict.keys():
                    for seg in sample_orgin_dict[chrom]:
                        sample_orgin_result.append(f'{seg["cytoband"]}({seg["origin"]})')
            # 过滤掉非嵌合 _nomos_ 保留含有嵌合的key
            sample_anno_info = {k: v for k, v in sample_anno_info_temp.items() if re.findall(r"_mos_", k)}
            x = ""
            each_roh = self.temp_roh.get(each_sample,"")
            if len(sample_anno_info) == 2:
                # print(sample_anno_info.keys())
                for i in sample_anno_info:
                    if '_no_report_gender' in i:
                        # 获取不带性别的指定分辨率下的注释信息
                        x = i
                logger.debug(f"{each_sample}->{x}")
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
                    if sample_orgin_result:
                        res3["CNV检测结果（来源）"] = ",".join(sample_orgin_result)
                    else:
                        res3["CNV检测结果（来源）"] = res3["CNV检测结果"]
                    res3["性别"] = res3.pop("Sex_chromosome_karyotype")
                    res3["结果解释"] = res3.pop("Description_CN")
                    # res3["结果说明"] = self.Format_结果说明(res3.pop("result_description"))
                    if not each_roh:
                        logger.debug(f"深测序roh为空, 直接导出浅测序结果注释")
                        res3["结果说明"] = self.Format_结果说明(res3.pop("result_description"))
                    else:
                        logger.debug(f"深测序roh不为空, 添加深测序roh注释进结果")
                        result_description_list = res3.pop("result_description")
                        new_result_description_list = self.add_snp_roh(result_description_list,each_roh)
                        res3["结果说明"] =self.Format_结果说明(new_result_description_list)
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
                    if sample_orgin_result:
                        res1["CNV检测结果（来源）"] = ",".join(sample_orgin_result)
                    else:
                        res1["CNV检测结果（来源）"] = res1["CNV检测结果"] 
                    res1["性别"] = res1.pop("Sex_chromosome_karyotype")
                    res1["结果解释"] = res1.pop("Description_CN")
                    # res1["结果说明"] = self.Format_结果说明(res1.pop("result_description"))
                    if not each_roh:
                        logger.debug(f"深测序roh为空, 直接导出浅测序结果注释")
                        res1["结果说明"] = self.Format_结果说明(res1.pop("result_description"))
                    else:
                        logger.debug(f"深测序roh不为空, 添加深测序roh注释进结果")
                        result_description_list = res3.pop("result_description")
                        new_result_description_list = self.add_snp_roh(result_description_list, each_roh)
                        res1["结果说明"] =self.Format_结果说明(new_result_description_list)
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
                logger.debug(f"{each_sample}->{x}")
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
                    if sample_orgin_result:
                        res3["CNV检测结果（来源）"] = ",".join(sample_orgin_result)
                    else:
                        res3["CNV检测结果（来源）"] = res3["CNV检测结果"] 
                    res3["性别"] = res3.pop("Sex_chromosome_karyotype")
                    res3["结果解释"] = res3.pop("Description_CN")
                    # res3["结果说明"] = self.Format_结果说明(res3.pop("result_description"))
                    if not each_roh:
                        logger.debug(f"深测序roh为空, 直接导出浅测序结果注释")
                        res3["结果说明"] = self.Format_结果说明(res3.pop("result_description"))
                    else:
                        logger.debug(f"深测序roh不为空, 添加深测序roh注释进结果")
                        result_description_list = res3.pop("result_description")
                        new_result_description_list = self.add_snp_roh(result_description_list, each_roh)
                        res3["结果说明"] =self.Format_结果说明(new_result_description_list)
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


def DoSomething_PGTM(zip_file_path, SNP, CNV, CNV_QC_,
                     platform_type, TCD, TCDIF, temp_roh,
                     NICS, NICS_data,
                     bp, ad,orgin_CNV_info):
    global logger
    project_name = Path(zip_file_path).name
    save_file1 = F"{zip_file_path}/Project_{project_name}_QC_Report.json"
    save_file2 = F"{zip_file_path}/Project_{project_name}_result_Report.json"
    if re.findall("_06AL_", project_name):
        logger.info(F"Code _06AL_ was recognized!")
        flag_chimeDetectionResult = True
    else:
        flag_chimeDetectionResult = False
    if platform_type == "NGS-Inst_Marsala":
        if not SNP and CNV:
            logger.info("仅含有CNV的质控数据")
            analysis_type = "CNV"
            if analysis_type == "CNV":
                with open(save_file1, 'w', encoding="utf-8") as F_CNV:
                    json.dump(CNV_QC_, F_CNV, sort_keys=True, ensure_ascii=False, indent="\t")
                if NICS:
                    # 不可能执行的部分
                    with open(save_file2, 'w', encoding="utf-8") as R_CNV:
                        json.dump(CNV_NewPipeReport(ad, TCD, TCDIF, temp_roh, NICS=NICS_data, flag_CDR=flag_chimeDetectionResult).Add_Breakpoint2SampleDict(), R_CNV, sort_keys=True, ensure_ascii=False, indent="\t")
                    return [save_file1, save_file2]
                else:
                    if re.findall("_MaReCs1|_PGTSR", project_name):
                        with open(save_file2, 'w', encoding="utf-8") as R_CNV:
                            json.dump(CNV_NewPipeReport(ad, TCD, TCDIF, temp_roh, NICS=None, bp=bp, flag_CDR=flag_chimeDetectionResult).Add_Breakpoint2SampleDict(), R_CNV, sort_keys=True, ensure_ascii=False, indent="\t")
                    elif re.findall("_PGTAH",project_name):
                        with open(save_file2, 'w', encoding="utf-8") as R_CNV:
                            json.dump(CNV_NewPipeReport(ad, TCD, TCDIF, temp_roh, NICS=None, flag_CDR=flag_chimeDetectionResult,orgin_dict=orgin_CNV_info).Add_Breakpoint2SampleDict(), R_CNV, sort_keys=True, ensure_ascii=False, indent="\t")
                    else:
                     
                        with open(save_file2, 'w', encoding="utf-8") as R_CNV:
                            json.dump(CNV_NewPipeReport(ad, TCD, TCDIF, temp_roh, NICS=None, flag_CDR=flag_chimeDetectionResult).Add_Breakpoint2SampleDict(), R_CNV, sort_keys=True, ensure_ascii=False, indent="\t")
                    return [save_file1, save_file2]
            else:
                logger.error("-a 参数值和包内容不符，脚本自动退出")
                sys.exit(1)
        else:
            logger.info("不包含SNP & CNV质控数据，稍后退出不作处理！")
            return []


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

# 读取orgin文件并解析为字典
def origin_to_json(file_path):
    result = {}
    print(file_path)
    if not os.path.exists(file_path):
        return {}
    
    # 染色体排序优先级：X, Y, 1-22
    CHROM_ORDER = {
        'X': 0, 'Y': 1, 
        **{str(i): i+1 for i in range(1, 23)}  # 1→2, 2→3, ..., 22→23
    }
    
    with open(file_path, 'r') as f:
        # 跳过注释行和标题行
        lines = [line.strip() for line in f if not line.startswith('#') and line.strip()]
        header = lines[0].split('\t')
        data_lines = lines[1:]

        for line in data_lines:
            cols = line.split('\t')
            embryo = cols[1]       # 第2列: Embryo (key1)
            chrom = cols[2]        # 第3列: Chromosome
            cnv = cols[5]          # 第6列: CNV (key2)
            conclusions = cols[12]  # 第13列: Conclusion (value)
            trans_conclusion = conclusions.split(":")
            transformed = []
            for conc in trans_conclusion:
                conc_str = conc.lower().replace("maternal","mat").replace("paternal","pat").replace("meiosisii","meiosis").replace("meiosisi","meiosis").replace("bph","meiosis").replace("sph","mitosis*")
                priority = 1 if conc_str.startswith(("mat", "pat")) else 0
                transformed.append((priority, conc_str, conc))
            # 排序：先其他(0)，后mat/pat(1)
            transformed.sort(key=lambda x: x[0])
            conclusions = ",".join([item[1] for item in transformed])
            

            if embryo not in result:
                result[embryo] = {}
            
            # 按染色体分组，并生成 cytoband 和 origin 结构
            if chrom not in result[embryo]:
                result[embryo][chrom] = []
            
            # 提取 cytoband（如 "-18" 或 "del(9)(p24.3p13.1)"）
            result[embryo][chrom].append({
                "cytoband": cnv,
                "origin": conclusions
            })
    
    # 对每个 Embryo 的染色体按指定顺序排序
    for embryo in result:
        sorted_chroms = sorted(
            result[embryo].keys(),
            key=lambda x: CHROM_ORDER.get(x, 100)  # 未知染色体排在最后
        )
        result[embryo] = {chrom: result[embryo][chrom] for chrom in sorted_chroms}
    
    return result

def Run_CNV_report(platform_type, CNV_Path, project_name, ZIP_res_dir, temp_roh,CNV_orgin_txt):
    # PGD & IB_PGD
    PGD_project = ['19AFSLASR', '19ASLASR', '19AFMSLASR', 'PGTMFSPFREE', 'PGD', 'PGTM', 'PGTMSR', 'PGTMFFREE',
                   'PGTMFREE', 'PGDRD2000', 'PGTMRD2000', 'PGDM', 'PGTMM', 'PGTMMRD2000',
                   'NIPGTM', 'PGDC', 'PGTMC', 'PGD_JX', 'HLA-PGD', 'HLAPGTM', 'HLA-PGDF',
                   'HLAPGTMF', 'PGTMCRD2000', 'PGDF', 'PGDFRD2000', 'PGTMFRD2000', 'PGTMF', 'PGTMFT',
                   'PGDF-S', 'PGTMFS', 'PGDFSP', 'PGDSR', 'PGTMFSP', 'PGDA', 'PGTMA', 'PBDSR',
                   'IB-PGD', 'IB_PGD', 'IBPGTM', 'IB']
    n = 0
    snp, cnv, NICS, bp = "", "", "", ""
    TCD, TCDIF, NICS_data = "", "", ""  # TCD=total_CNV_detail， TCDIF=total_CNV_detail_IF
    SNP_QC_res, CNV_QC_res, MT_dd_list, json_BP_CNV = {}, {}, {}, {}
    json_QC_CNV, json_CNV_info = {}, {}
    orgin_CNV_info = origin_to_json(CNV_orgin_txt)
    withXY_mos_1M = {}
    noXY_mos_1M = {}
    Chromosome_Diagram = defaultdict(list)
    logger.debug(F"开始遍历目录: {CNV_Path}")
    
    def traverse_directory(directory_path, file_paths_):
        path = Path(directory_path)
        for item in path.iterdir():
            if item.is_file():
                file_paths_.append(item)  # 将文件路径添加到列表中
            elif item.is_dir():
                traverse_directory(item, file_paths_)  # 递归遍历子文件夹
    
    file_paths = []
    traverse_directory(CNV_Path, file_paths)
    sub_files = file_paths
    
    logger.debug(len(sub_files))
    
    for each in sub_files:
        if F"{project_name}_total_data_sts.json" in each.name:
            cnv = each
            logger.info("检查到包含CNV的质控数据!")
            with open(each) as f_QC_CNV:
                json_QC_CNV = json.loads(f_QC_CNV.read())
            n += 1
        elif F"{project_name}_total_cnv.json" in each.name:
            with open(each) as f_CNV:
                json_CNV_info = json.loads(f_CNV.read())
            n += 1
        elif "_withXY_mos.png" in each.name:
            key_name = Path(each).name.replace("_withXY_mos.png", "")
            withXY_mos_1M[key_name] = each.as_posix().replace(ZIP_res_dir, "")
            n += 1
        elif "_noXY_mos.png" in each.name:
            key_name = Path(each).name.replace("_noXY_mos.png", "")
            noXY_mos_1M[key_name] = each.as_posix().replace(ZIP_res_dir, "")
            n += 1
        elif "_cytoband.png" in each.name:
            key_name = Path(each).name.replace("_cytoband.png", "")
            Chromosome_Diagram[key_name].append(each.as_posix().replace(ZIP_res_dir, ""))
        elif "info.xlsx" in each.name and contains(each.as_posix(), "local_report_Source"):
            # TODO PGD&IB_PGD相关项目需要单独处理，其info表格只有3列，其他都是5列
            #  -- 2022年4月28日 -- v9_bin/generate_info_and_fullgraph.py
            # print(each)
            pro_code = project_name.split("_")[1]
            if pro_code in PGD_project:
                info_sample = pd.read_excel(each, header=None, na_values="-", dtype=str,
                                            na_filter=False, keep_default_na=False,
                                            names=["样本条码", "样本名称", "CNV结果"],
                                            sheet_name='info', index_col="样本条码").to_dict(orient="index")
            else:
                info_sample = pd.read_excel(each, header=None, na_values="-", dtype=str,
                                            na_filter=False, keep_default_na=False,
                                            names=["样本条码", "样本名称", "性别", "CNV结果", "女方"],
                                            sheet_name='info', index_col="样本条码").to_dict(orient="index")
            n += 1
        # TODO 注释读取_total_cnv.json和意外发现注释都直接读_Additional_cnv_details.tsv
        elif "_Additional_cnv_details.tsv" in each.name:
            TCDIF = pd.read_csv(each, header=0, delimiter="\t")
            n += 1
        elif F"{project_name}_data.tsv" in each.name:
            n += 1
            MT_df = pd.read_csv(each, header=0, delimiter="\t", index_col="##SampleID", usecols=["##SampleID", "MT_CN"])
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
        # elif "_breakpoint_by_project.json" in each:
        elif "_BreakPoints_Kmeans.json" in each.name:
            if re.findall("_MaReCs|_PGTSR", project_name):
                with open(each) as BP_CNV:
                    json_BP_CNV = json.loads(BP_CNV.read())
                pro, bp = parse_breakpoint_info(json_BP_CNV)
            else:
                bp = None
            n += 1
    
    if cnv:
        additional_data = CNV_SupplementData(info_sample, withXY_mos_1M, noXY_mos_1M, Chromosome_Diagram)
    else:
        additional_data = {}
    if platform_type == "NGS-Inst_Marsala":
        CNV_QC_res = CNV_QC(json_QC_CNV, withXY_mos_1M, json_CNV_info, MT_dd_list).ReadNGS_CVN_QC2NewNGS_QC_Report()
    else:
        logger.error("当前脚本只处理NGS的Marsala-Inst_pgta数据!")
        sys.exit(1)
    
    if n > 0:
        # TODO Do Something
        TCD_json_CNV_info = json_CNV_info.copy()
        QC_file = DoSomething_PGTM(ZIP_res_dir, snp, cnv, CNV_QC_res,
                                   platform_type, TCD_json_CNV_info, TCDIF, temp_roh,
                                   NICS, NICS_data,
                                   bp, additional_data,orgin_CNV_info)
        test_info = '\n'.join(QC_file)
        logger.info(F"qc4LIMS json files: \n{test_info}")
    else:
        logger.warning("检查输入文件没有检查到CNV的结果文件，n为0!!!")
        sys.exit(1)


def arguments():
    parser = argparse.ArgumentParser(prog=__describe__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''
                                     \033[0;31;33m
                                     该脚本用于NGS-Inst_Marsala
                                     PGTA|MaReCs|PGTSR|PGTAH等不涉及连锁项目的QC&results文件整理
                                     \033[0m\n'''))
    parser.add_argument('--version', action='version', version=__describe__)
    parser.add_argument('--ResultPath', '-r', dest='ZP', type=str, default=None, required=True,
                        help='Results file path, '
                             'e.g. /data/Project_2023Q3/Project_YKSZ_PGTAH_230715_21F_01_0836/pgtah_analysis/result '
                             '\n /data/Project_2023Q3/Project_YKSZ_PGTAH_230715_21F_01_0836')
    parser.add_argument('--PlatformType', '-p', dest='PT', type=str, default="NGS-Inst_Marsala", choices=['NGS-Inst_Marsala'],
                        help='NGS-Inst_Marsala 分析逻辑方可使用当前脚本进行处理!')
    # if len(sys.argv) == 1:
    #     parser.print_help(sys.stderr)
    #     sys.exit(1)
    return parser.parse_args()


def  get_roh_info(each_file):
    if not Path(each_file).exists():
        logger.error(f"{each_file} 不存在， 退出")
        sys.exit(1)
    logger.debug(each_file)
    try:
        df_temp = pd.read_excel(each_file, header=1, index_col="Sample", dtype=str,na_values="", na_filter=False, keep_default_na=False)
    except Exception as ERROR:
        logger.error(ERROR)
        df_temp = pd.read_excel(each_file, header=0, index_col="Sample", dtype=str,na_values="", na_filter=False, keep_default_na=False)
    info_temp = df_temp.to_dict(orient='index')
    roh_dd = {}
    for k, v in info_temp.items():
        roh_anno = v.get("ROHanno", "")
        if roh_anno:
            roh_dd[k] = roh_anno
    return roh_dd


def main():
    global logger
    
    args = arguments()
    input_results_ = args.ZP  # /data/Project_2023Q3/Project_YKSZ_PGTAH_230715_21F_01_0836/pgtah_analysis/result
    platform_type = args.PT  # NGS-Inst_Marsala
    CNV_zip = []
    if "_analysis/result" in input_results_:
        input_results = input_results_
        cnv_run_dir = Path(input_results_).parent / "CNV"
    else:
        # 分析目录或者结果目录都可以
        input_results = F"{Path(input_results_).as_posix()}/pgtah_analysis/result"
        cnv_run_dir =  Path(input_results_) / "pgtah_analysis"/ "CNV"
    
  
    
    for each in Path(input_results).glob("**/CNV/*.zip"):
        logger.debug(F"CNV: {each}")
        CNV_zip.append(each.as_posix())
    
    if CNV_zip:
        input_zip = CNV_zip[0]
        CNV_Path = Path(input_zip).parent.as_posix()
        # /data/Project_2023Q3/Project_YKSZ_PGTAH_230715_21F_01_0836/pgtah_analysis/result/YKSZ_PGTAH_230715_21F_01_0836/CNV/
        ZIP_res_dir = Path(CNV_Path).parents[0].as_posix()
        # /data/Project_2023Q3/Project_YKSZ_PGTAH_230715_21F_01_0836/pgtah_analysis/result/YKSZ_PGTAH_230715_21F_01_0836
        # 1) 解压CNV包,删除压缩包
        extract_zipfile(input_zip, CNV_Path)
        if Path(input_zip).exists():
            Path(input_zip).unlink()
    else:
        logger.error(F"没有检查到CNV/*.zip文件,请检查是否分析出CNV包!")
        temp_local = []
        for each in Path(input_results).glob("**/CNV/local_report_Source"):
            logger.debug(F"CNV-local_report_Source: {each}")
            temp_local.append(each.as_posix())
        if not temp_local:
            logger.error(F"没有检查到CNV/local_report_Source目录,退出!")
            sys.exit(1)
        
        CNV_Path = Path(temp_local[0]).parent.as_posix()
        ZIP_res_dir = Path(CNV_Path).parents[0].as_posix()
    
    version_json = F"{input_results}/version.json"
    if Path(version_json).exists():
        res_ver = shutil.copy(version_json, ZIP_res_dir)
        logger.info(F"版本信息文件： {res_ver}")
    else:
        logger.info(F"无版本信息JSON")
    
    # /data/Project_2023Q3/Project_YKSZ_PGTAH_230715_21F_01_0836/pgtah_analysis/result/YKSZ_PGTAH_230715_21F_01_0836/CNV/YKSZ_PGTAH_230715_21F_01_0836.zip
    project_name = Path(CNV_Path).parents[0].name
    CNV_orgin_txt = f"{input_results_}/pgtah_analysis/{project_name}/CNVorigin/{project_name}_origin.txt"
    tag_ld_pgta = cnv_run_dir / project_name / "ld_pgta_analysis" / "ld_pgta_run_time_stats.json"
    snp2cnv = f"{input_results}/{project_name}/{project_name}.snp2cnv.cn.xlsx"
    temp_roh = get_roh_info(snp2cnv)
    logger.debug(f"{temp_roh=}")

    if not tag_ld_pgta.exists():
        logger.error(f"ld-pgta可能跑断, 请检查后重试！ 没有 {tag_ld_pgta=}")
        sys.exit(1)
    if re.findall("_06AL_", project_name):
        cmd_add = f"/data/home/zhouyajun/miniconda3/bin/python  /data/biomed/sPGD_Project_Temp/Script_to_be_online/Individualized_Demand_Scripts/06AL_CNV_SNP_data_v1.py   -d  {input_results_} "
        logger.debug(f"执行: {cmd_add}")
        ExecuteCMD(cmd_add)
    
    if re.findall("_Control_", project_name, flags=re.IGNORECASE):
        logger.error("该流程Control项目直接不处理！")
        sys.exit(0)
    elif re.findall("PGTA|MaReCs", project_name):
        logger.info(f"识别到项目：{project_name}")
        # TODO CNV 部分功能输出,并且输出JSON后解压压缩包到CNV目录
        # 2) 调整后的CNV路径进行直接读取整合
        Run_CNV_report(platform_type, CNV_Path, project_name, ZIP_res_dir,temp_roh,CNV_orgin_txt)
        save_file = F"Project_{project_name}_QC_Report.json"
        QC_file = F"{ZIP_res_dir}/{save_file}"
        if Path(QC_file).exists():
            logger.info(F"QC_Report.json已经输出准备导入ASA_QC_Data")
            # 3) SNP部分逻辑需要追加QC部分
            PGTAH_Type(ZIP_res_dir).DoSomethingPGTAH()
            # 4) 重压缩所有文件, zip包存放至 Path(ZIP_res_dir).parents[2].as_posix()
            zip_file = F"{Path(ZIP_res_dir).parents[2].as_posix()}/Project_{project_name}.zip"
            ZIP_result(ZIP_res_dir, zip_file)
        else:
            logger.error(F"没有检查到QC_Report.json, 退出!")
            sys.exit(1)
    
    else:
        logger.error(F"项目代码不识别!")
        sys.exit(1)


def extract_zipfile(zipfile_path, destination_folder):
    with zipfile.ZipFile(zipfile_path, 'r', compression=zipfile.ZIP_DEFLATED) as target:
        target.extractall(destination_folder)
    logger.info("CNV/*.zip 解压完毕!")


def ZIP_result(zip_dir, zip_package_file):
    import os
    
    def add_folder_to_zip(zip_file_, folder_path, parent_folder=''):
        for item in os.listdir(folder_path):
            item_path = os.path.join(folder_path, item)
            if os.path.isfile(item_path):
                zip_file_.write(item_path, arcname=os.path.join(parent_folder, item))
            elif os.path.isdir(item_path):
                add_folder_to_zip(zip_file_, item_path, os.path.join(parent_folder, item))
    
    with zipfile.ZipFile(zip_package_file, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zip_file:
        if Path(zip_dir).exists():
            try:
                add_folder_to_zip(zip_file, zip_dir)
            except Exception as e:
                logger.error(e)


class PGTAH_Type(object):
    """
    PGTAH的信息储存在XX_result_Report.json
    """
    
    def __init__(self, ZIP_res_dir):
        self.ZIP_res_dir = ZIP_res_dir
        self.pro_name = Path(self.ZIP_res_dir).name
    
    @staticmethod
    def remove_suffix(input_string, suffix):
        # removesuffix只有python3.9才有，向下兼容构建当前函数
        if input_string.endswith(suffix):
            return input_string[:-len(suffix)]
        return input_string
    
    def Get_snp2cnv(self):
        CSV_file = F"{self.ZIP_res_dir}/{self.pro_name}.snp2cnv.cn.xlsx"
        if not Path(CSV_file).exists():
            logger.error(F"文件未检查到: {CSV_file}")
            sys.exit(1)
        
        logger.info(F"Get: {CSV_file}")
        # df = pd.read_csv(CSV_file, comment="#", header=0, usecols=["sample", "result"], index_col="sample")
        try:
            df = pd.read_excel(CSV_file, header=1, index_col="Sample", dtype=str,
                               na_values="", na_filter=False, keep_default_na=False)
        except Exception as ERROR:
            logger.error(ERROR)
            df = pd.read_excel(CSV_file, header=0, index_col="Sample", dtype=str,
                               na_values="", na_filter=False, keep_default_na=False)
        info_dd = df.to_dict(orient='index')
        new_dd = {}
        for k, v in info_dd.items():
            new_dd[k] = v['Heteroploid_result']
        return new_dd, info_dd
    
    def DoSomethingPGTAH(self, res_file="QC_Report"):
        dd = defaultdict(dict)
        all_sample = []
        uneven_dir = F"{self.ZIP_res_dir}/uneven"
        if Path(uneven_dir).exists():
            for each in Path(uneven_dir).iterdir():
                if each.name.endswith(".vaf.density_peak.png"):
                    # 后续异常染色体会单独放入uneven,优先尝试匹配2位染色体，不行就1位
                    if re.findall(r"(?:chr[1-9][0-9]|chr[1-9|X|Y])", Path(each).name):
                        logger.debug(F"跳过该样本记录(包含染色体异常图)：{each}")
                        continue
                    # 默认全部样本都有
                    all_sample.append(self.remove_suffix(Path(each).name, ".vaf.density_peak.png"))
                    dd["BAF_png"][self.remove_suffix(Path(each).name, ".vaf.density_peak.png")] = [each.as_posix().replace(Path(self.ZIP_res_dir).as_posix(), "")]
                elif each.name.endswith("_LOHScore.png"):
                    dd["LOH_PDF"][self.remove_suffix(Path(each).name, "_LOHScore.png")] = [each.as_posix().replace(Path(self.ZIP_res_dir).as_posix(), "")]
                elif each.name.endswith("_UnevenScore.png"):
                    dd["UN_PDF"][self.remove_suffix(Path(each).name, "_UnevenScore.png")] = [each.as_posix().replace(Path(self.ZIP_res_dir).as_posix(), "")]

        llr_pngs = F"{self.ZIP_res_dir}/CNV/"
        dd_pngs = {}
        for each_upd in Path(llr_pngs).glob("**/llr_pngs/*.chrom_llr.upd.png"):
            sp_name = Path(each_upd).name.replace(".chrom_llr.upd.png", "")
            dd_pngs[sp_name] = [each_upd.as_posix().replace(Path(self.ZIP_res_dir).as_posix(), "")]

        CNV_Thumbnail = defaultdict(list)
        for PNG_each in Path(llr_pngs).glob("**/other_check_graph/*.png"):
            if not PNG_each.name.endswith("_CNV.png"):
                barcode_name = PNG_each.name.split("_chr")[0]  # _chr分割默认
                CNV_Thumbnail[barcode_name].append(PNG_each.as_posix().replace(Path(self.ZIP_res_dir).as_posix(), ""))
        
        info_dd, all_dd = self.Get_snp2cnv()
        save_file = F"Project_{self.pro_name}_{res_file}.json"
        QC_file = F"{self.ZIP_res_dir}/{save_file}"
        with open(QC_file, "r", encoding="utf-8") as JF:
            dd_temp = json.load(JF)
        
        result = {}
        for key in all_sample:
            # logger.debug(F"sample_barcode: {key}")
            sample_dd = dd_temp[key]
            BAF = dd["BAF_png"].get(key, [])
            # 2023年9月26日 和南飞协调改成SNP_QC_Data 2023年9月26日改回ASA_QC_Data
            # 2023年12月7日 根据报告端需求调整输出图片内容,三体输出vaf+uneven。upd输出loh+uneven
            Heteroploid_result = info_dd.get(key, "")
            # TODO 2024年2月1日 -- 经过IT协调调整到_QC_Report.json，字段放入到ASA_QC_Data中，占用Key：VARIATION_RESULT_List: {{"pic_name": ""}}
            CNV_Thumbnail_pngs = CNV_Thumbnail.get(key, '')
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
                logger.error(F"检查到不被识别的类型: {key} --> {CNV_Thumbnail_pngs}")
                # sys.exit(1)
            
            if re.findall("三倍体|Triploid", Heteroploid_result, flags=re.IGNORECASE):
                BAF_value = BAF
                uneven_value = dd["UN_PDF"].get(key, [])
                loh_value = []
            elif re.findall("单亲二倍体|UPD", Heteroploid_result, flags=re.IGNORECASE):
                BAF_value = []
                uneven_value = dd["UN_PDF"].get(key, [])
                loh_value_temp = dd["LOH_PDF"].get(key, [])
                if not loh_value_temp:
                    # 搜索“.chrom_llr.upd.png”文件在CNV目录，递归检索
                    loh_value = dd_pngs.get(key, [])
                else:
                    loh_value = loh_value_temp
            else:
                # logger.info(F"{key} --> 非'三体' or UPD 不提供VAF or UPD点图")
                logger.info(F"{key} --> 非'三体' or UPD 不提供VAF or loh点图，强制全部PGTAH释放UnevenScore.png")
                # 2024年2月23日 增加 “PGTAH项目正常和异常的样本在json中都放出uneven这个图”
                BAF_value = []
                # uneven_value = []
                uneven_value = dd["UN_PDF"].get(key, [])
                loh_value = []

            t2 = {"BAF图谱": BAF_value, "亲本污染分析": Heteroploid_result,
                  "异倍体图谱(Uneven)": uneven_value,
                  "异倍体图谱(LOH)": loh_value,
                  "VARIATION_RESULT_List": dict(temp_CNV_info).get("VARIATION_RESULT_List",
                                                                   {"pic_name": "",
                                                                    "Mutation_type": "",
                                                                    "geneName": "",
                                                                    "variationInfo": ""})}
            logger.debug(F"测试：{key} | --> {t2}")
            # t2 = {"BAF图谱": BAF, "亲本污染分析": info_dd.get(key, ""),
            #       "异倍体图谱(Uneven)": dd["UN_PDF"].get(key, []),
            #       "异倍体图谱(LOH)": dd["LOH_PDF"].get(key, [])}
            t2.update(all_dd[key])
            comb_dd = self.MergeDict(sample_dd, {"ASA_QC_Data": t2})
            # sample_dd.update(t2)
            result[key] = comb_dd
        
        # 合并CNV结果再在输出
        with open(QC_file, "w", encoding="utf-8") as JF:
            json.dump(result, JF, sort_keys=True, ensure_ascii=False, indent="\t")
        logger.info(F"新数据覆盖完毕: {Path(QC_file).name}")
    
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


if __name__ == "__main__":
    sys_platform = platform.system()
    HostName = socket.gethostname()
    format_time = time.strftime('%Y_%m_%d|%H:%M:%S', time.localtime(time.time()))
    logger = create_logging('LIMS2ReportSystem_Inst-Marsala')
    start_info = f"Current Run Platform：{sys_platform}||Current Host：{HostName}"
    logger.info(start_info)
    logger.info("Start：{}".format(__describe__))
    # TODO Check Results
    # 测试路径: /data/home/zhouyajun/TestData/20230728_pgta_complete/Project_YKSZ_PGTAH_230715_21F_01_0836/
    main()
