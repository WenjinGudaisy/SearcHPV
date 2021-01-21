import os
import sys
import pysam
import pandas as pd
import pickle

# #read combined dic
# res = f'{out_dir}/call_fusion_virus/combinedContigDic.pickle'
# with open(res,'rb') as inputFile:
#     combinedContigDic = pickle.load(inputFile)
def read_mapping_info(out_dir):
    #read mapping info
    sitesPath = f'{out_dir}/assemble/'
    listSites = os.listdir(sitesPath)
    contigSamDictHPV = {}
    contigSamDictGenome = {}
    contigLengthDict = {}
    contigSRDict = {}
    if listSites != ['']:  
        for site in listSites:
            if ".sh" not in site:
                hpv_bam = f"{out_dir}/call_fusion_virus/{site}/{site}.contigToHPV.sort.bam"
                genome_bam = f"{out_dir}/call_fusion_virus/{site}/{site}.contigToGenome.sort.bam"
                try:
                    samFile = pysam.AlignmentFile(hpv_bam,'rb')
                except:
                    print("Contig cound not mapping to HPV",site)
                    continue
                    

                for read in samFile:
                    readList = []
                    contigName = read.qname
                    key = site+"."+contigName
                    if read.cigar != []:
                        readList.append(read.query_alignment_length) #matched BP
                        readList.append(read.cigarstring) #cigar
                        readList.append(read.pos) #pos
                        readList.append(read.cigar) #cigar turple

                        if key in contigSamDictHPV.keys():
                            contigSamDictHPV[key].append(readList)
                        else:
                            contigSamDictHPV[key] = [readList]
                            
                            
                try:
                    samFile = pysam.AlignmentFile(genome_bam,'rb')
                except:
                    print("Contig cound not mapping to Genome",sample,site)
                    continue
                    

                for read in samFile:
                    readList = []
                    contigName = read.qname
                    key = site+"."+contigName
                    
                    if read.cigar != []:
                        readList.append(read.query_alignment_length) #matched BP
                        readList.append(read.cigarstring) #cigar
                        readList.append(read.pos) #pos
                        readList.append(read.cigar) #cigar turple

                        if key in contigSamDictGenome.keys():
                            contigSamDictGenome[key].append(readList)
                            
                        else:
                            contigSamDictGenome[key] = [readList]

            #read length from contig file
            
                contigPath = f"{out_dir}/assemble/{site}/pearOutput/{site}.all.fa.cap.contigs"

                with open(contigPath) as contigFile:
                    contigs = contigFile.read().rstrip().split('>')
                    for each in contigs:
                        if each != '':
                            rows = each.rstrip().split('\n')
                            key = site+"."+rows[0]
                            allSeq = ''
                            for seq in rows[1:]:
                                allSeq += seq
                                #allSeq = allSeq.replace("N","")
                                contigLengthDict[key] = len(allSeq)#full length

                #read sr num from srNum.txt
                srPath = f'{out_dir}/call_fusion_virus/{site}/srNum.txt'
                with open(srPath) as sr:
                    info = sr.read().rstrip().split('\n')
                    for eachinfo in info:
                        infoList = eachinfo.split('\t')
                        key = infoList[0]
                #                     print(key)
                        try:
                            srNum = infoList[1]
                            contigSRDict[key] = int(srNum)
                        except:
                            print("didn't find srNUM.txt")
                            print(key)
                

    #print(contigLengthDict)    
    #combine information together
    contigDict = {}
    for key,value in contigSamDictHPV.items():
        if key in contigSamDictGenome.keys():
            try:
                contigDict[key] = [contigSamDictHPV[key],contigSamDictGenome[key],contigLengthDict[key],contigSRDict[key]]
            except:
                print(key,value)
                sys.exit("Lack of information.")
    return contigDict


#calculate hpv insertion
# find insertion point of HPV
#for regular contig hpv insertion point between h/s and m
#for regular contig, HPV insertion point between H and S
def cal_hpv_ins(contigDict,out_dir):
    sitesPath = f'{out_dir}/assemble/'
    listSites = os.listdir(sitesPath)
    contigRegular = {}

    i = 0
    for key,value in contigDict.items():
            #for those contigs has 1 HPV alignment hit
        if len(value[0]) == 1:
                contigRegular[key] = value
                i += 1
    for key,value in contigRegular.items():
        if value[0][0][3][0][0] == 0:#start with match, insertion point is pos + matched length
            hpvInsertion = value[0][0][2] + value[0][0][0]
        elif value[0][0][3][0][0] == 4 or 5:#start with clip, insertion point is pos
            hpvInsertion = value[0][0][2]
        value.append(hpvInsertion)

    #for circle contig
    contigCircle = {}

    i = 0
    for key,value in contigDict.items():
            #for those contigs has 1 HPV alignment hit
        if len(value[0]) > 1:
            if value[0][0][2] == 0:#pos equal to 0
                #print(key,value)
                cigarTurpleDic = {}
                matchLength = int(value[0][0][0]) + int(value[0][1][0])#need add to new dic
                clipStart = value[0][0][-1][0][1]
                for each in value[0][1][-1]:
                    if each[0] == 0:
                        matchEnd = each[1]
                if clipStart <= matchEnd or value[0][0][-1][0] == 5:
                    newCigarTurple = []
                    for each in value[0][1][-1]:
                        newCigarTurple.append(list(each))
                    i = 0
                    for each in value[0][1][-1]:
                        if each[0]== 0:
                            newCigarTurple[i][1] = matchLength
                            break
                        i += 1
                    print(value,newCigarTurple)
                    # if newCigarTurple[i+1][0] == 4 or newCigarTurple[i+1][0] == 5:
                    #     newCigarTurple[i+1][1] = newCigarTurple[i+1][1] - int(value[0][0][0])
                    # else:
                    #     print("Error",newCigarTurple)

                    newValue = [[[matchLength]+value[0][1][1:-1] + [newCigarTurple]],value[1],value[-2],value[-1]]
                else:
                    newCigarTurple = []
                    for each in value[0][0][-1]:
                        newCigarTurple.append(list(each))
                    newStart = clipStart - matchEnd
                    #print(newCigarTurple)
                    i = 0
                    for each in value[0][0][-1]:
                        if each[0]== 0:
                            newCigarTurple[i][1] = matchLength
                        i += 1
                    if newCigarTurple[0][0] == 4 or newCigarTurple[0][0] == 5:
                        #print('1',newCigarTurple[0])
                        newCigarTurple[0][1] = newCigarTurple[0][1] - int(value[0][1][0])
                    else:
                        print("Error",newCigarTurple)
                    newValue = [[[matchLength]+[value[0][0][1]] + [newStart] + [newCigarTurple]],value[1],value[-2],value[-1]]
                contigCircle[key] = newValue
            #print(newValue)
            

#             for eachHit in value[0]:#hit mapping to hpv
#                 for each in eachHit[-1]:#cigar turple: change cigar turple as well for choose contig
#                     if each[1] in cigarTurpleDic.keys():
#                         cigarTurpleDic[each[1]].append(each[0]) #let the cigar number be key, cigar type be value
#                     else:
#                         cigarTurpleDic[each[1]] = [each[0]]
#             print(cigarTurpleDic)
#             for key1,value1 in cigarTurpleDic.items():
#                 if len(value1) >= 2:#if cigar number are same, the lengh is 2
#                     overLap = key1
            
#             #find an element in turple with overlap as cigar number and H or S as cigar type
#             listNeed = []
#             for eachHit in value[0]:
#                 listTurple = []
#                 for each in eachHit[-1]:
#                     listEach = list(each)
#                     listTurple.append(listEach)
#                 for each in listTurple:
#                     if len(listTurple) == 3:
#                         if each[1] == overLap and (each[0] == 4 or each[0] == 5):#change it as match
#                             each[0] = 0
#                             listNeed = listTurple
# #                     print(listNeed)
#             if len(value[0][0][3]) > len(value[0][1][3]):
#                 newValue = [[[matchLength]+value[0][0][1:-1] + [listNeed]],value[1],value[-2],value[-1]]
#             else:
#                 newValue = [[[matchLength]+value[0][1][1:-1] + [listNeed]],value[1],value[-2],value[-1]]
#             contigCircle[key] = newValue# only leave the first hit, all the information is in txt file

#calculate hpv insertion
# find insertion point of HPV
#for circle contig hpv insertion point between h/s and m 
#but there are 2 slot of m
    for key,value in contigCircle.items():
        if value[0][0][3][0][0] == 0:#start with match, insertion point is 0+second matched piece
            hpvInsertion = 0 + value[0][0][3][1][1]
        elif value[0][0][3][0][0] == 4 or 5:#start with clip, insertion point is pos
            hpvInsertion = value[0][0][2]
        value.append(hpvInsertion)
            

#for wierd contigs
    contigWierd = {}
    for site in listSites:
        if ".sh" not in site:
            for key,value in contigDict.items():
        #         if key == 'Sample_109356.11.38317802.Contig1':
                if site in key:
                    #3 situation, 
                        #1: due to HPV circle structure, pos is 0 and number of matched bp is equal to the mate's
                        #soft/hard clip length
                        #2: due to dupilication of HPV. matched region has overlap region with mate's matched region
                        #3: matched region has no overlap region, seperated with mate's matched region

                        #for 1, add the length of this read and its mate together, combine two list to be one

                    #2,3 wierd hits need check
                    if len(value[0]) > 1:
                        if value[0][0][2] == 0:
                            pass
                        else:
                            contigWierd[key] = value



    #calculate hpv insertion
    # find insertion point of HPV
    #find matched region from genome mapping, find a most close softclip region in hpv hit
    #print(contigWierd)
    for key,value in contigWierd.items():
    #     if key == 'Sample_109356.11.38317802.Contig1':
        if len(value[1]) == 1:
            matchedGenome = value[1][0][0]
            diff = 7000#a big number
            
            for each in value[0]:
                for slot in each[3]:
                    #print(slot)
                    if int(slot[0]) == 4 or int(slot[0]) == 5:
                        clipHPV = slot[1]
                        #print(clipHPV)
                        if diff > abs(clipHPV-matchedGenome):
                            diff = abs(clipHPV-matchedGenome)
                            #print(diff)
                            wantSlot = each
                            targetClip = clipHPV
                
            #print('1',wantSlot,targetClip,matchedGenome)
            i = 0
            if wantSlot[3][0][1] == targetClip:
                hpvInsertion = wantSlot[2]#pos
            elif wantSlot[3][-1][1] == targetClip:
                hpvInsertion = wantSlot[2] + wantSlot[0]#pos + matchedLength
            #else:
                #print('2',wantSlot,targetClip,matchedGenome)
    #             value = value + [hpvInsertion]

            value.append(hpvInsertion)
        else:

            pos = int(key.split('.')[1])
            diff = pos
            valueWant = value[1][0]
            for each in value[1]:
                if diff > abs(each[2]-pos):
                    valueWant = each
            #print(valueWant)
            matchedGenome = valueWant[0]
            for each in value[0]:
                for slot in each[3]:
                    if int(slot[0]) == 4 or int(slot[0]) == 5:
                        clipHPV = slot[1]
                        if diff > abs(clipHPV-matchedGenome):
                            diff = abs(clipHPV-matchedGenome)
                            wantSlot = each
                            targetClip = clipHPV
            #print('3',wantSlot,targetClip,matchedGenome)
            i = 0
            if wantSlot[3][0][1] == targetClip:
                hpvInsertion = wantSlot[2]#pos
            elif wantSlot[3][-1][1] == targetClip:
                hpvInsertion = wantSlot[2] + wantSlot[0]#pos + matchedLength
            #else:
                #print('4',wantSlot,targetClip,matchedGenome)
            #print("wantSlot",wantSlot)
    #             value[2] = [wantSlot]
            #print(hpvInsertion)
            #print(value)
    #             value = value + [hpvInsertion]
            #print(value)
            value.append(hpvInsertion)

    #combine three apart together
    combinedContigDic = {}
    for key,value in contigRegular.items():
        combinedContigDic[key] = value
    for key,value in contigCircle.items():
        combinedContigDic[key] = value
    for key,value in contigWierd.items():
        combinedContigDic[key] = value
    return combinedContigDic



def select_contig(combinedContigDic):
    newSiteList = []
    for each in combinedContigDic.keys():
        newSite = '.'.join(each.split('.')[:-1])
        newSiteList.append(newSite)

    newSiteList = set(newSiteList)
    #print(newSiteList)
    #select for contigs
    selectedAllContig = {}
    contigStatus = {}

    for site in newSiteList:
        maxLengthLeft = 0
        maxSRLeft = 0
        maxLengthRight = 0
        maxSRRight = 0
        numContig = 0
        for key,value in combinedContigDic.items():
            if site in key:
                numContig += 1
        if numContig == 1:# if only one contig in site
            for key,value in combinedContigDic.items():
                if site in key:
                    matchedbp = 0
                    #for contigs has mutiple alignments to HPV
                    for each in value[0]:
                        matchedbp += each[0]
                    if float(matchedbp)/float(value[-3]) >0.1 and float(matchedbp)/float(value[-3]) <0.95 and value[-2] >= 10: #SR > 10 and matchedBP/Length > 0.3 < 1
                        selectedAllContig[site] = [[key.split('.')[-1]] + value]
                        contigStatus[site] = [[key.split('.')[-1]],'one']
    #                     print('one',key,value)
                    else:
                        selectedAllContig[site] = [[key.split('.')[-1]] + value + ['lowConfidence']]
                        contigStatus[site] = [[key.split('.')[-1]],'one low']
        #                     print('one low',key,value)
        else:# if there are more than 1 contig
            #count max lenght left and right for high confidence reads
            for key,value in combinedContigDic.items():
                if site in key:
                    matchedbp = 0
                    for each in value[0]:
                        matchedbp += each[0]
                    isLeft = False
                    
                    for i,each in enumerate(value[1][0][3][:-1]):
                        
                        if each[0] == 0 and (value[1][0][3][i+1][0] == 4 or value[1][0][3][i+1][0] == 5):
                            
                            isLeft = True
                    if isLeft and float(matchedbp)/float(value[-3]) > 0.1 and float(matchedbp)/float(value[-3]) < 0.95 and value[-2] >= 10 : # if cigar mapping to genome start with match count max length and mac sr
                        if int(value[-3]) > maxLengthLeft and float(matchedbp)/float(value[-3]) != 1:# trim out the contig fully matched to hpv
                            maxLengthLeft = int(value[-3])
                            
                        
                        # except:
                        #     print("Error",key,value)
                    elif not isLeft and float(matchedbp)/float(value[-3]) > 0.1 and float(matchedbp)/float(value[-3]) < 0.95 and value[-2] >= 10: # if cigar mapping to genome start with clip count max
                        if int(value[-3]) > maxLengthRight and float(matchedbp)/float(value[-3]) != 1:
                            maxLengthRight = int(value[-3])
                            
                            
                            #print(value[-3])

                        # except:
                        #     print("Error",key,value)
            
                

            #count max sr num left and right, do not need to identify confidence because these two parameters will only used when no high confidence contigs exist
            for key,value in combinedContigDic.items():
                if site in key:
                    matchedbp = 0
                    for each in value[0]:
                        matchedbp += each[0]
                    isLeft = False
                    
                    for i,each in enumerate(value[1][0][3][:-1]):
                        
                        if each[0] == 0 and (value[1][0][3][i+1][0] == 4 or value[1][0][3][i+1][0] == 5):
                            
                            isLeft = True
                    if isLeft: # if cigar mapping to genome start with match count max length and mac sr
                        

                        if int(value[-2]) > maxSRLeft:

                            maxSRLeft = int(value[-2])
                            
                        # except:
                        #     print("Error",key,value)
                    else: # if cigar mapping to genome start with clip count max
                        
                        if int(value[-2]) > maxSRRight:
                            maxSRRight = int(value[-2])

            #print(site,maxLengthLeft,maxLengthRight,maxSRLeft,maxSRRight)
            #select contigs for high confidence contigs
            addLeft = False
            addRight = False
            
            for key,value in combinedContigDic.items():
                if site in key:
                    matchedbp = 0
                    for each in value[0]:
                        matchedbp += each[0]
                    isLeft = False
                    
                    for i,each in enumerate(value[1][0][3][:-1]):
                        
                        if each[0] == 0 and (value[1][0][3][i+1][0] == 4 or value[1][0][3][i+1][0] == 5):
                            #print(each)
                            isLeft = True
                    
                    if  int(value[-3]) == maxLengthLeft and float(matchedbp)/float(value[-3]) > 0.1 and float(matchedbp)/float(value[-3]) < 0.95 and isLeft and value[-2] >= 10:
                        #print("left",key,value)
                        
                        if site in selectedAllContig.keys():
                            selectedAllContig[site].append([key.split('.')[-1]] + value)
                            contigStatus[site].append([[key.split('.')[-1]],'left high'])
                            addLeft = True
                            break
                        else:
                            selectedAllContig[site] = [[key.split('.')[-1]] + value]
                            contigStatus[site] = [[key.split('.')[-1]],'left high']
                            addLeft = True
                            break
            for key,value in combinedContigDic.items():
                if site in key:
                    matchedbp = 0
                    for each in value[0]:
                        matchedbp += each[0]
                    isLeft = False
                    
                    for i,each in enumerate(value[1][0][3][:-1]):
                        
                        if each[0] == 0 and (value[1][0][3][i+1][0] == 4 or value[1][0][3][i+1][0] == 5):
                            #print(each)
                            isLeft = True
                    if  int(value[-3]) == maxLengthRight and float(matchedbp)/float(value[-3]) > 0.1 and float(matchedbp)/float(value[-3]) < 0.95 and not isLeft and value[-2] >= 10:
                        
                        #print("right",key,value)
                        if site in selectedAllContig.keys():
                            selectedAllContig[site].append([key.split('.')[-1]] + value)
                            contigStatus[site].append([[key.split('.')[-1]],'right high'])
                            addRight = True
                            break
                        else:
                            selectedAllContig[site] = [[key.split('.')[-1]] + value]
                            contigStatus[site] = [[key.split('.')[-1]],'right high']
                            addRight = True
                            break
            #add contigs from left for low confidence contigs
            
            if addLeft == False:
                for key,value in combinedContigDic.items():
                    if site in key:
                        matchedbp = 0
                        for each in value[0]:
                            matchedbp += each[0]
                        isLeft = False

                        for i,each in enumerate(value[1][0][3][:-1]):
                            if each[0] == 0 and (value[1][0][3][i+1][0] == 4 or value[1][0][3][i+1][0] == 5):
                                #print(each)
                                isLeft = True
                        if isLeft and int(value[-2]) == maxSRLeft:
                            #print("left",key,value)
                            if site in selectedAllContig.keys():
                                selectedAllContig[site].append([key.split('.')[-1]] + value + ['lowConfidence'])
                                contigStatus[site].append([[key.split('.')[-1]],'left low'])
                                addLeft = True
                                break
                            else:
                                selectedAllContig[site] = [[key.split('.')[-1]] + value + ['lowConfidence']]
                                contigStatus[site] = [[key.split('.')[-1]],'left low']
                                addLeft = True
                                break

            #add contigs from right for low confidence contigs
            if addRight == False:
                for key,value in combinedContigDic.items():
                    if site in key:
                        matchedbp = 0
                        for each in value[0]:
                            matchedbp += each[0]
                        isLeft = False
                        
                        for i,each in enumerate(value[1][0][3][:-1]):
                            if each[0] == 0 and (value[1][0][3][i+1][0] == 4 or value[1][0][3][i+1][0] == 5):
                                #print(each)
                                isLeft = True
                        if not isLeft and int(value[-2]) == maxSRRight:

                            #print("right",key,value)
                            if site in selectedAllContig.keys():
                                selectedAllContig[site].append([key.split('.')[-1]] + value + ['lowConfidence'])
                                contigStatus[site].append([[key.split('.')[-1]],'right low'])
                                addRight = True
                                break
                            else:
                                selectedAllContig[site] = [[key.split('.')[-1]] + value + ['lowConfidence']]
                                contigStatus[site] = [[key.split('.')[-1]],'right low']
                                addRight = True
                                break

            if addLeft != True and addRight != True:
                print("error",site,value)
    return selectedAllContig
            

    #     j = 0
    #     #if could not find best contigs use criteria above
    #     if i == numContig:
    #         for key,value in combinedContigDic.items():
    #             if site in key:
    #                 matchedbp = 0
    #                 for each in value[0]:
    #                     matchedbp += each[0]
    #                 #select high quality contig first
    #                 if site in selectedAllContig.keys():
    #                     if float(matchedbp)/float(value[-3]) >0.1 and float(matchedbp)/float(value[-3]) <0.95 and value[-2] >= 10:
    #                         selectedAllContig[site].append([key.split('.')[-1]] + value)
    #                         contigStatus[site].append([[key.split('.')[-1]],'maxonly in high'])
    #                         #print(site,[key.split('.')[-1]] + value)
    #                     else:
    #                         j += 1
    # #                             selectedAllContig[site].append([key.split('.')[-1]] + value + ['low Confidence'])
    #                 else:
    #                     if float(matchedbp)/float(value[-3]) >0.1 and float(matchedbp)/float(value[-3]) <0.95 and value[-2] >= 10:
    #                         #print(site,[key.split('.')[-1]] + value)
    #                         selectedAllContig[site]= [[key.split('.')[-1]] + value]
    #                         contigStatus[site] = [[key.split('.')[-1]],'maxonly in high']
    #                     else:
    #                         j += 1
    # #                             selectedAllContig[site] = [[key.split('.')[-1]] + value + ['low Confidence']]
    #     #if no high quality contigs, select low quality contigs with maximum supportive reads
    #     if j == numContig:
    #         maxSR = 0
    #         for key,value in combinedContigDic.items():
    #             #print(key,site)
    #             if site in key:
    #                 if value[-2] > maxSR: 
    #                     maxSR = value[-2]
    #                     wantContig = [key.split('.')[-1]] + value + ['lowConfidence']
    #                     #print("low",wantContig)
    # #                             print(key,value)
    
    #         if site in selectedAllContig.keys():
    #             selectedAllContig[site].append(wantContig)
    #             contigStatus[site].append([[key.split('.')[-1]],'low'])
    #         else:
    #             selectedAllContig[site]= [wantContig]
    #             contigStatus[site] = [[key.split('.')[-1]],'low']

# for each in combinedContigDic:
#     if 'Sample_109371.X.71608238' in each:
#         print(each,combinedContigDic[each])
# print(selectedAllContig['Sample_109371.X.71608238'])

###############
#extract siteConfidence
#out_dir: output directory for seacHPV
def siteConf(out_dir):
    siteConfidence = {}
    with open(f'{out_dir}/call_fusion/all.filtered.clustered.result') as site_res:
        siteList = site_res.read().split(';')
        #print(siteList)
        for eachSite in siteList:
            chro = eachSite.split(':')[0]
            site =  eachSite.split(':')[1]
            siteConfidence[chro + '.' + site] = [eachSite.split(':')[2],eachSite.split(':')[3],eachSite.split(':')[-1]]
    return siteConfidence

#############
#extract contig sequence
#out_dir:
def extractContigSeq(contigDic,out_dir):
    contigSeq = {}
    for key,value in contigDic.items():
        chro = key.split('.')[0]
        if 'GL' not in chro:
            file = f'{out_dir}/assemble/{key}/pearOutput/{key}.all.fa.cap.contigs'
            with open(file) as contigs:
                contigs = contigs.read().split('>')
        #         print(contigs)

                for contig in contigs:
                    if contig != '':
                        contig = contig.rstrip('')
                        contigName = contig.split('\n')[0]
                        for selectedContig in value:
                            selectedContigName = selectedContig[0]
                            if selectedContigName == contigName:
                                sequence = ''.join(contig.split('\n')[1:])
                                if key in contigSeq:
                                    contigSeq[key].append([contigName,sequence])
                                else:
                                    contigSeq[key] = [[contigName,sequence]]
    return contigSeq
        
##############
# check if repetive region in contig
def check_rep(contigSeq,key):
    repFlag = False
    repA = 'A'*80
    repT = 'T'*80
    repC = 'C'*80
    repG = 'G'*80
    for eachSeq in contigSeq[key]:
        if repA in eachSeq[1] or repT in eachSeq[1] or repC in eachSeq[1] or repG in eachSeq[1]:
            repFlag = True
    return repFlag


################
#filter results with only one low confident contig and low number of supportive reads, repetitive region
def filter_res(selectedAllContig,siteConfidence,out_dir):
    contigSeq = extractContigSeq(selectedAllContig,out_dir)
    filteredSelectedContig = {}
    
    for key in selectedAllContig:
        if int(siteConfidence[key][0]) < 4 and int(siteConfidence[key][1]) < 4:
            i = 0
            for contig in selectedAllContig[key]:
                i += 1
        #if two contigs
            if i > 1:
                filteredSelectedContig[key] = selectedAllContig[key]
            else:
        #only one contig
                #if low confident contig
                if contig[-1] == 'lowConfidence':
                    #if supprotive read greater than 3
                    if contig[-3] > 3:
                        if check_rep(contigSeq,key)== False:
                            filteredSelectedContig[key] = selectedAllContig[key]
                else:
                    if check_rep(contigSeq,key) == False:
                        filteredSelectedContig[key] = selectedAllContig[key]
        else:
            if check_rep(contigSeq,key) == False:
                filteredSelectedContig[key] = selectedAllContig[key]
    return filteredSelectedContig



################
#write output to file
#out_dir: output directory for seacHPV
def write_to_file(filteredSelectedContig,siteConfidence,out_dir):
    outputPath = f'{out_dir}/call_fusion_virus/'
    if os.path.isfile(outputPath + 'HPVfusionPointContig.txt'):
        os.system(f'rm {outputPath}/HPVfusionPointContig.txt')
    
            
    with open(outputPath + 'HPVfusionPointContig.txt','w') as output: 
        output.write('contigName\tgenomeInsertionChromosome\tgenomeInsertionPoint\thpvInsertionPoint\tSiteConfidence\tContigConfidence\n')
    #     output.write('sampleName\tcontigName\tgenomeInsertionChromosome\tgenomeInsertionPoint\thpvInsertionPoint\tConfidence\n')
        for key,value in filteredSelectedContig.items():
    #         if key != 'Sample_81279.3.189612850':
            
            chro = key.split('.')[0]
            if 'GL' not in chro:
                genomeInsertion = key.split('.')[1]
                siteConf = siteConfidence[chro + '.' + genomeInsertion][2]
                if siteConf == 'high':
                    siteConfOutput = 'High Confidence Site'
                else:
                    siteConfOutput = 'Low Confidence Site'
                for eachContig in value:
                    #print(key,eachContig)
                    contigName = key + '.'+eachContig[0]
                    if eachContig[-1] != 'lowConfidence':
                        hpvInsertion = eachContig[-1]
                        contigConfidence = 'High Confidence Contig'
                    else:
                        hpvInsertion = eachContig[-2]
                        contigConfidence = 'Low Confidence Contig'
                    output.write(f'{contigName}\t{chro}\t{genomeInsertion}\t {hpvInsertion}\t{siteConfOutput}\t{contigConfidence}\n')
                output.write('\n')


    #extract sequence for contigs and write to file
    outputFile = outputPath + "ContigsSequence.fa"
    with open(outputFile,'w') as output:
        for key,value in filteredSelectedContig.items():
            chro = key.split('.')[0]
            if 'GL' not in chro:
                file = f'{out_dir}/assemble/{key}/pearOutput/{key}.all.fa.cap.contigs'
                with open(file) as contigs:
                    contigs = contigs.read().split('>')
            #         print(contigs)

                    for contig in contigs:
                        if contig != '':
                            contig = contig.rstrip('')
                            contigName = contig.split('\n')[0]
                            for selectedContig in value:
                                selectedContigName = selectedContig[0]
                                if selectedContigName == contigName:
                                    sequence = ''.join(contig.split('\n')[1:])
                                    output.write('>' + key + '.' + contigName + '\n')
                                    output.write(sequence + '\n')
                    output.write('\n')



    f = open(outputPath +  'filteredSelectedContig.pickle',"wb")
    pickle.dump(filteredSelectedContig,f)
    f.close()

    # f = open(outputPath +  'weirdContig.pickle',"wb")
    # pickle.dump(contigWierd,f)
    # f.close()

    # f = open(outputPath + 'combinedContigDic.pickle','wb')
    # pickle.dump(combinedContigDic,f)
    # f.close()

    return None

