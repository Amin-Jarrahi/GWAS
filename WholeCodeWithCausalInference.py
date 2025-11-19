
# # Libraries we will use

from tabulate import tabulate
import pandas as pd
import numpy as np
import random
from scipy import stats
from sklearn.ensemble import RandomForestRegressor
import sklearn
from sklearn.feature_selection import RFE
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.neighbors import KNeighborsRegressor
from sklearn.linear_model import LinearRegression
import time
import matplotlib.pyplot as plt
import dowhy
import dowhy.plotter
from dowhy import CausalModel
import os

def BuildVCFandHeight(Genename):
#VCF DATA
    filepath = 'D:/users/ajarrah1/Desktop/HeightGenes/AllGenes/'+Genename+'.vcf.gz'
    vcfdata = pd.read_csv(filepath,compression='gzip', header=3385, sep='\t')
    fullvcfdata = vcfdata
    vcfdata=vcfdata.replace(to_replace ="0|0", value = 0.0)
    vcfdata=vcfdata.replace(to_replace ="0|1", value = 1.0)
    vcfdata=vcfdata.replace(to_replace ="1|0", value = 1.0)
    vcfdata=vcfdata.replace(to_replace ="1|1", value = 2.0)
    vcfdata=vcfdata.replace(to_replace =".|.", value = 3.0)
    columnnames = vcfdata.columns.to_list()
    vcfdata.columns = [columnnames]
 #Height DATA
    df = pd.read_excel('D:/Users/ajarrah1/Desktop/HeightGenes/GTEx-Height-Demographics.xlsx')
    df=df.replace(to_replace ="20-29", value = 0.0)
    df=df.replace(to_replace ="30-39", value = 0.0)
    df=df.replace(to_replace ="40-49", value = 1.0)
    df=df.replace(to_replace ="50-59", value = 2.0)
    df=df.replace(to_replace ="60-69", value = 3.0)
    df=df.replace(to_replace ="70-79", value = 4.0)
    df=df.replace(to_replace ="1", value = 2.0)
    df=df.replace(to_replace ="2", value = 1.0)

    
    print(df)
    XlsxToDic = df.to_dict()
    dic = {}
    sex_dic = {}
    age_dic = {}
    #print(len(a["SUBJID"]))
    for i in range(len(XlsxToDic["SUBJID"])):
        key = XlsxToDic['SUBJID'][i]
        value = XlsxToDic['HGHT'][i]
        dic.update({key: value})
    for sex in range(len(XlsxToDic["SUBJID"])):
        key_s = XlsxToDic['SUBJID'][sex]
        value_s = XlsxToDic['SEX'][sex]
        sex_dic.update({key_s: value_s})
    for age in range(len(XlsxToDic["SUBJID"])):
        key_a = XlsxToDic['SUBJID'][age]
        value_a = XlsxToDic['AGE'][age]
        age_dic.update({key_a: value_a})
    #print(dic)
    #print(len(dic))
    targetdatalist = list()
    sexdata = list()
    agedata = list()
    for key1 in vcfdata.head():
        key1 = str(key1)
        key1 = ''.join(str(key1).split(','))
        key1 = ''.join(str(key1).split(')'))
        key1 = ''.join(str(key1).split('('))
        key1 = ''.join(str(key1).split('"'))
        key1 = ''.join(str(key1).split("'"))
    
        if key1 in dic.keys():
            targetdatalist.append(float(dic[key1]))
            sexdata.append(float(sex_dic[key1]))
            agedata.append(float(age_dic[key1]))
    
    currenttarget = pd.DataFrame(targetdatalist, columns=['Height'])
    currenttarget["Sex"] = sexdata
    currenttarget["Age"] = agedata
    print(currenttarget)

    vcfdata=vcfdata.drop(['#CHROM', 'POS', 'ID','REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], axis=1)
    vcfdata = vcfdata.transpose()
    print(vcfdata.shape)
    return vcfdata, currenttarget

def BuildVCF(Genename):
#VCF DATA
    filepath = 'D:/users/ajarrah1/Desktop/HeightGenes/AllGenes/'+Genename+'.vcf.gz'
    vcfdata = pd.read_csv(filepath,compression='gzip', header=3385, sep='\t')
    fullvcfdata = vcfdata
    vcfdata=vcfdata.replace(to_replace ="0|0", value = 0.0)
    vcfdata=vcfdata.replace(to_replace ="0|1", value = 1.0)
    vcfdata=vcfdata.replace(to_replace ="1|0", value = 1.0)
    vcfdata=vcfdata.replace(to_replace ="1|1", value = 2.0)
    vcfdata=vcfdata.replace(to_replace =".|.", value = 3.0)
    columnnames = vcfdata.columns.to_list()
    vcfdata.columns = [columnnames]
    vcfdata=vcfdata.drop(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], axis=1)
    vcfdata = vcfdata.transpose()
    print(vcfdata.shape)
    return vcfdata


def BuildSack(GeneData, n, sack=None):
#Concatinating the DATA
    if n==0:
        sack = GeneData  
    else:
        sack = pd.concat([sack, GeneData], axis=1)
    print(sack.shape)
    return sack

def BuildNDArrayOFGeneDataAndHeightData(sack,HeightData):
    #sack=sack.drop(["ID"], axis=1)
    data = sack.values
    data = data[:, :].astype(str)
    data=data
    X = data
    #HeightData.drop(['Sex', 'Age'], axis=0)
    print(HeightData.shape)
    Y = HeightData['Height'].values.astype('float32')
    Y = Y.ravel()
    X = np.asarray(X).astype('float32')
    print(X.shape,Y.shape)
    return X, Y

def BuildNDArrayOFSack(filtered_sack):
    #sack=sack.drop(["ID"], axis=1)
    print("filtered_sack", filtered_sack.shape)
    data = filtered_sack.values
    data = data[:, :].astype(str)
    data=data
    X = data
    X = np.asarray(X).astype('float32')
    print("BuildNDArrayOFSack X.shape", X.shape)
    np.savetxt('D:/Users/ajarrah1/Desktop/HeightGenes/WholeCodeWithCausalInferenceResults/errorsss'+str(numberoffeatures)+'.csv',X,delimiter=",")
    return X

#pd.set_option('display.max_columns', None)
#pd.set_option('display.max_rows', None)
#pd.set_option('display.max_colwidth', None)
#genesack += vcfdata
#print(genesack)
    

def CVGen(Indexnum, n):  
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(Indexnum), n):  
        yield Indexnum[i:i + n] #yielding the groups of 140 members from the shuffled list


# ## Define ML methods

def RFFeatureSelection(TrainVariants,TrainHeight):
    regr = RandomForestRegressor(n_estimators =100, random_state=42)
    regr.fit(TrainVariants, TrainHeight)
    importances = regr.feature_importances_
    return importances


def RFand2SRF(TrainVariants,TrainHeight,TestVariants,TestHeight, WholeTrainVariants,WholeTrainHeight):
    reg =RandomForestRegressor(n_estimators =100, random_state=42).fit(TrainVariants,TrainHeight)
    y_test_pred=reg.predict(TestVariants)

    reg1 = RandomForestRegressor(n_estimators =100, random_state=42).fit(WholeTrainHeight.reshape(-1,1), WholeTrainVariants)
    x_2sls_pred = reg1.predict(y_test_pred.reshape(-1,1))

    reg2 = RandomForestRegressor(n_estimators =100, random_state=42).fit(x_2sls_pred, TestHeight)
    y_test_pred_causal = reg2.predict(x_2sls_pred)
    return y_test_pred,y_test_pred_causal

def RFReg(TrainVariants,TrainHeight,TestVariants,TestHeight,numcom):
    regr = RandomForestRegressor(n_estimators =100, random_state=42)
    regr.fit(TrainVariants, TrainHeight)
    PredictedHeight = regr.predict(TestVariants)
    importances = regr.feature_importances_
    return PredictedHeight, importances

def KNNand2SKNN(TrainVariants,TrainHeight,TestVariants,TestHeight, WholeTrainVariants,WholeTrainHeight):
    reg = KNeighborsRegressor(n_neighbors = 5).fit(TrainVariants,TrainHeight)
    y_test_pred=reg.predict(TestVariants)

    reg1 = KNeighborsRegressor(n_neighbors = 5).fit(WholeTrainHeight.reshape(-1,1), WholeTrainVariants)
    x_2sls_pred = reg1.predict(y_test_pred.reshape(-1,1))

    reg2 = KNeighborsRegressor(n_neighbors = 5).fit(x_2sls_pred, TestHeight)
    y_test_pred_causal = reg2.predict(x_2sls_pred)
    return y_test_pred,y_test_pred_causal

def OLSand2SLS(TrainVariants,TrainHeight,TestVariants,TestHeight, WholeTrainVariants,WholeTrainHeight):
    reg =LinearRegression().fit(TrainVariants,TrainHeight)
    y_test_pred=reg.predict(TestVariants)

    reg1 = LinearRegression().fit(WholeTrainHeight.reshape(-1,1), WholeTrainVariants)
    x_2sls_pred = reg1.predict(y_test_pred.reshape(-1,1))
    reg2 = LinearRegression().fit(x_2sls_pred, TestHeight)
    y_test_pred_causal = reg2.predict(x_2sls_pred)
    return y_test_pred,y_test_pred_causal

def Causal_Model(X,HeightData, numberoffeatures):
    snp = np.arange(0,len(X.T))
    snp = snp.tolist()
    snp = [str(i) for i in snp]
    df_of_X = pd.DataFrame(X, columns = snp)
    df_of_X.reset_index()
    HeightData.reset_index()
    df_of_X_merged = pd.concat([df_of_X,HeightData], axis=1)
    #df_of_X_merged=df_of_X_merged.drop(['Age','Sex'], axis=1)
   
    #print(df_of_X_merged.head())
   # print(df_of_X_merged.shape)
    CIColumnsList = ["ID",
                    "Causal Effect Estimate",
                    "New Effect Random Common Cause",
                    "Random Common Cause P-Value",
                    "New Effect Random placebo",
                    "Random placebo P-Value",
                    "New effect of removing random sub set of data",
                    "removing random sub set of data P-Value" ]


    CI = pd.DataFrame(columns=CIColumnsList)

    for treatment in snp:
        model = CausalModel(
            data=df_of_X_merged,
            treatment=treatment, #causal_model
            outcome='Height',
            common_causes=['Age','Sex'],#common causes are variables that affects both treatment and outcome. These variable are also caleed confounders 
            intruments=None #it affects outcome only with association with exposure of interest(treatment)
            )
        
            #identify causal effectusing properties of the formal causal gragh
        identified_estimand = model.identify_effect(proceed_when_unidentifiable=True)
        #print("identified_estimand: ",identified_estimand)

        #estimate the causal effect
        estimate = model.estimate_effect(
                identified_estimand,
                method_name='backdoor.linear_regression',
                test_significance=True
                )
        #print("estimate: ", estimate)
        


        #adding a random common cause variable
        res_random=model.refute_estimate(identified_estimand,
                                        estimate,
                                        method_name='random_common_cause'
                                        )
        #print("res_random: ", res_random)
        file_sack = open('txtfolder/file_sack.txt','w+')
        file_sack.write(str(res_random))
        file_sack.close()
        with open('txtfolder/file_sack.txt','r') as file_sack:
            lines = file_sack.readlines()

        x=['New effect', 'p value']
        res_random_results=[]
        for line in lines:
            parts = line.split(":")
            if parts[0]in x:
                res_random_results.append(round(float(parts[1].split("\\")[0]),7))
            
        #print(res_random_results)
        file_sack.close()
        os.remove('txtfolder/file_sack.txt')

        #replacing treatment with random (placebo) variable
        res_placebo=model.refute_estimate(identified_estimand,
                                        estimate,
                                        method_name='placebo_treatment_refuter',
                                        placebo_type='permute'
                                        )
       # print("res_placebo: ",res_placebo)
        file_sack = open('txtfolder/file_sack.txt','w+')
        file_sack.write(str(res_placebo))
        file_sack.close()
        with open('txtfolder/file_sack.txt','r') as file_sack:
            lines = file_sack.readlines()
            
        x=['New effect', 'p value']
        res_placebo_results=[]
        substrings=[]
        for line in lines:
            parts = line.split(":")
            if parts[0]in x:
                res_placebo_results.append(round(float(parts[1].split("\\")[0]),7))

                   
       # print(res_placebo_results)
        file_sack.close()
        os.remove('txtfolder/file_sack.txt')
        
        #removing random sub set of data
        res_subset= model.refute_estimate(identified_estimand,
                                        estimate,
                                        method_name='data_subset_refuter',
                                        subset_fraction=0.9
                                        )
        #print("res_subset: ", res_subset)
        file_sack = open('txtfolder/file_sack.txt','w+')
        file_sack.write(str(res_subset))
        file_sack.close()
        with open('txtfolder/file_sack.txt','r') as file_sack:
            lines = file_sack.readlines()
            
        x=['New effect', 'p value']
        res_subset_results=[]
        substrings=[]
        for line in lines:
            parts = line.split(":")
            if parts[0]in x:
                res_subset_results.append(round(float(parts[1].split("\\")[0]),7))
            
       # print(res_subset_results)
        file_sack.close()
        os.remove('txtfolder/file_sack.txt')

        CI_Results = {"ID": treatment,
                    "Causal Effect Estimate": round(estimate.value,7),
                    "New Effect Random Common Cause": res_random_results[0],
                    "Random Common Cause P-Value": res_random_results[1],
                    "New Effect Random placebo": res_placebo_results[0],
                    "Random placebo P-Value": res_placebo_results[1],
                    "New effect of removing random sub set of data": res_subset_results[0],
                    "removing random sub set of data P-Value": res_subset_results[1]
                      }

        CI = pd.concat([CI, pd.DataFrame(data=CI_Results, index=[0])], ignore_index=True)
#        print('CI Results', CI_Results)
 #       print(CI)
    CI.to_csv('D:/Users/ajarrah1/Desktop/HeightGenes/WholeCodeWithCausalInferenceResults/CI '+str(numberoffeatures)+'.csv')
    sorted_CI=CI.sort_values(by='Causal Effect Estimate', ascending=False)

    sorted_CI.to_csv('D:/Users/ajarrah1/Desktop/HeightGenes/WholeCodeWithCausalInferenceResults/CISorted '+str(numberoffeatures)+'.csv')

    
            



# ## Define Other Helpful Things

def FilterData(imp, numberoffeatures):
    RFImportanceData = imp
    print("RFImportanceData", RFImportanceData.shape)

    #RFImportanceData_filtered = RFImportanceData[RFImportanceData[0]>=cutoffvalue]
    
    RFImportanceData.columns = ["Importances"]
    RFImportanceData.index.name = "ID"
    sorted_RFImportanceData=RFImportanceData.sort_values(by='Importances', ascending=False)
    RFImportanceData_filtered = sorted_RFImportanceData[:numberoffeatures]
    print("RFImportanceData_filtered", RFImportanceData_filtered.shape)

    return RFImportanceData_filtered


## #  dropping data based on RFImportances

def DropUnimportantData(RFImportanceData_filtered, sack):

    columnnames = sack.columns.to_list()
    print("sack", sack)
    print("sack", sack.shape)

    RFImportanceData_filtered = RFImportanceData_filtered.T
    RFImp_columnName = RFImportanceData_filtered.columns.to_list()
    print('RFImp_columnName', len(RFImp_columnName))
    RFImp_columnNameList=[]
    
    for RFImp_col_name in RFImp_columnName:
        RFImp_col_name = str(RFImp_col_name)
        RFImp_col_name = ''.join(str(RFImp_col_name).split(','))
        RFImp_col_name = ''.join(str(RFImp_col_name).split(')'))
        RFImp_col_name = ''.join(str(RFImp_col_name).split('('))
        RFImp_col_name = ''.join(str(RFImp_col_name).split('"'))
        RFImp_col_name = ''.join(str(RFImp_col_name).split("'"))
        RFImp_columnNameList.append(RFImp_col_name)
    print('RFImp_columnNameList', len(RFImp_columnNameList))

    drop_list = list()
    keep_list = list()
    for col_name in columnnames:
        col_name = str(col_name)
        col_name = ''.join(str(col_name).split(','))
        col_name = ''.join(str(col_name).split(')'))
        col_name = ''.join(str(col_name).split('('))
        col_name = ''.join(str(col_name).split('"'))
        col_name = ''.join(str(col_name).split("'"))
        if col_name not in RFImp_columnNameList:
            drop_list.append(col_name)
        if col_name in RFImp_columnNameList:
            keep_list.append(col_name)
            
    print('drop_list', drop_list)
    print('keep_list', keep_list)

    sack = sack[keep_list]

 #   sack=sack.drop(drop_list, axis=1)
    print(sack.shape)
    return sack

def cvgroups(X,Y,CV):
    TV = np.delete(X, CV, axis=0)
    TE = np.delete(Y, CV, axis=0)
    TeV = X[CV]
    TeE = Y[CV]
    return TV,TE,TeV,TeE

def getR2(TE,Pred,OldR2):
    tss = sum(TE**2)
    rss = sum((TE - Pred)**2)
    r2 = 1 - rss/tss
    OldR2.append(r2)
    return OldR2

def getMAE(TE,Pred,OldMAE):
    mae=mean_absolute_error(TE, Pred)
    OldMAE.append(mae)
    return OldMAE

def getMSE(TE,Pred,OldMSE):
    mse=mean_squared_error(TE, Pred)
    OldMSE.append(mse)
    return OldMSE

def IDNameList(Genename, IDList):
    filepath = 'D:/users/ajarrah1/Desktop/HeightGenes/AllGenes/'+Genename+'.vcf.gz'
    IDdata = pd.read_csv(filepath,compression='gzip', header=3385, sep='\t', usecols=["ID"])
    
    IDdata = IDdata.transpose()
    IDdata = IDdata.to_dict()
    for keys in IDdata:
        IDList.append(IDdata[keys]["ID"])
    return IDList



random.seed(2022)
Indexnum = np.arange(0,838,1) #producing a list of numbers from 0 to 838(number of people)
random.shuffle(Indexnum)       #Shuffling the list
CVGroups = list(CVGen(Indexnum, 140))   #deviding the shuffled list into groups of 140 members



DFColumnsList = ['Number of Features',
                     'OLS_R2','STD_OLS_R2','2SLS_R2','STD_2SLS_R2',
                     'RF_R2','STD_RF_R2','2SRF_R2','STD_2SRF_R2',
                     'KNN_R2','STD_KNN_R2','2SKNN_R2','STD_2SKNN_R2',
                     'OLS_MAE','STD_OLS_MAE','2SLS_MAE','STD_2SLS_MAE',
                     'RF_MAE','STD_RF_MAE','2SRF_MAE','STD_2SRF_MAE',
                     'KNN_MAE','STD_KNN_MAE','2SKNN_MAE','STD_2SKNN_MAE',
                     'OLS_MSE','STD_OLS_MSE','2SLS_MSE','STD_2SLS_MSE',
                     'RF_MSE','STD_RF_MSE','2SRF_MSE','STD_2SRF_MSE',
                     'KNN_MSE','STD_KNN_MSE','2SKNN_MSE','STD_2SKNN_MSE']

df = pd.DataFrame(columns=DFColumnsList)


# ## Run ML Models

    ##Imprtance DataFrame  

df2 = pd.DataFrame(columns=DFColumnsList)

#Genelist = ['ENSG00000143476.17', 'ENSG00000146197.8','ENSG00000065600.12',
 #               'ENSG00000119681.11',  'ENSG00000149257.13', 'ENSG00000157766.15',
  #              'ENSG00000159899.14', 'ENSG00000169047.5','ENSG00000182752.9',
   #             'ENSG00000185960.13','ENSG00000116288.12','ENSG00000145335.15',
    #        'ENSG00000177628.15','ENSG00000185345.18','ENSG00000188906.14']
Genelist = ['ENSG00000143476.17', 'ENSG00000146197.8',
                'ENSG00000119681.11',  'ENSG00000149257.13', 'ENSG00000157766.15',
                'ENSG00000159899.14', 'ENSG00000169047.5','ENSG00000182752.9',
                'ENSG00000185960.13']
#Genelist = ['ENSG00000143476.17']
print(len(Genelist))
n=0
IDList = []
ID_List = []
for gene in Genelist:

    if n==0:
        GeneData, HeightData = BuildVCFandHeight(gene)
    else:
        GeneData = BuildVCF(gene)
    ID_List = IDNameList(gene, IDList)
    ID_Set = set(ID_List)
    print("ID_List ",len(ID_List))
    print("ID_Set ",len(ID_Set))
    if n==0:
        sack = BuildSack(GeneData, n)
    else:
        sack = BuildSack(GeneData, n, sack)
    n += 1
        
impColumnsList = [ID_List]
imp = pd.DataFrame(columns=impColumnsList)

#sack.to_pickle('D:/Users/ajarrah1/Desktop/HeightGenes/WholeResults/ConcatinatedSNPs.pkl')
X,Y = BuildNDArrayOFGeneDataAndHeightData(sack,HeightData)
sack.columns=impColumnsList
importances = RFFeatureSelection(X,Y)

imp.loc[len(imp.index)] = importances
imp = imp.T

imp.to_csv('D:/Users/ajarrah1/Desktop/HeightGenes/WholeCodeWithCausalInferenceResults/impOutput.csv')

    
#Cutoffvaluerange=np.linspace(0.0001, 0.002, 20)
#reversed_cutoffvaluerange=np.flip(Cutoffvaluerange)
#print(reversed_cutoffvaluerange)
numberoffeatures_list = [2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,
                         200,300,400,500,600,700,800,900,1000,1500,2000,2500,
                         3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,
                         8500,9000,9500,10000,12000,14000,16000,18000,20000]
for numberoffeatures in numberoffeatures_list:
        

    RFImportanceData_filtered = FilterData(imp,numberoffeatures)
    filtered_sack = DropUnimportantData(RFImportanceData_filtered, sack)
    X = BuildNDArrayOFSack(filtered_sack)
    Causal_Model(X,HeightData,numberoffeatures)




    R2ValuesRF = []
    MAEValuesRF = []
    MSEValuesRF = []
    R2Values2SRF = []
    MAEValues2SRF = []
    MSEValues2SRF = []
    R2ValuesOLS = []
    MAEValuesOLS = []
    MSEValuesOLS = []
    R2Values2SLS = []
    MAEValues2SLS = []
    MSEValues2SLS = []
    R2ValuesKNN = []
    MAEValuesKNN = []
    MSEValuesKNN = []
    R2Values2SKNN = []
    MAEValues2SKNN = []
    MSEValues2SKNN = []
    print("MSEValues2SKNN", len(MSEValues2SKNN))

    WholeTrainVariants = X
    print("X", X.shape)
    WholeTrainHeight = Y.ravel()
    print("Y", Y.shape)
    
    for CV in CVGroups:
            #deviding into 4 groups of TrainVariants,TrainHeight,TestVariants and TestHeight
        TrainVariants,TrainHeight,TestVariants,TestHeight = cvgroups(X,Y.ravel(),CV)
            #Predict the test set RF and 2SRF
        y_test_pred_RF,y_test_pred_causal_2SRF = RFand2SRF(TrainVariants,TrainHeight,
                                                           TestVariants,TestHeight,
                                                           WholeTrainVariants,WholeTrainHeight)
            #Calculate the R2      
        R2ValuesRF = getR2(TestHeight,y_test_pred_RF,R2ValuesRF)
            #Calculate the MAE      
        MAEValuesRF = getMAE(TestHeight,y_test_pred_RF,MAEValuesRF)
            #Calculate the MSE      
        MSEValuesRF = getMSE(TestHeight,y_test_pred_RF,MSEValuesRF)
            #Calculate the R2 Causal     
        R2Values2SRF = getR2(TestHeight,y_test_pred_causal_2SRF,R2Values2SRF)
            #Calculate the MAE Causal     
        MAEValues2SRF = getMAE(TestHeight,y_test_pred_causal_2SRF,MAEValues2SRF)
            #Calculate the MSE Causal     
        MSEValues2SRF = getMSE(TestHeight,y_test_pred_causal_2SRF,MSEValues2SRF)
            #Predict the test set KNN and 2SKNN
        y_test_pred_KNN,y_test_pred_causal_2SKNN = KNNand2SKNN(TrainVariants,TrainHeight,
                                                               TestVariants,TestHeight,
                                                               WholeTrainVariants,WholeTrainHeight)
            #Calculate the R2      
        R2ValuesKNN = getR2(TestHeight,y_test_pred_KNN,R2ValuesKNN)
            #Calculate the MAE      
        MAEValuesKNN = getMAE(TestHeight,y_test_pred_KNN,MAEValuesKNN)
            #Calculate the MSE      
        MSEValuesKNN = getMSE(TestHeight,y_test_pred_KNN,MSEValuesKNN)
            #Calculate the R2 Causal     
        R2Values2SKNN = getR2(TestHeight,y_test_pred_causal_2SKNN,R2Values2SKNN)
            #Calculate the MAE Causal     
        MAEValues2SKNN = getMAE(TestHeight,y_test_pred_causal_2SKNN,MAEValues2SKNN)
            #Calculate the MSE Causal     
        MSEValues2SKNN = getMSE(TestHeight,y_test_pred_causal_2SKNN,MSEValues2SKNN)
            #Predict the test set OLS and 2SLS
        y_test_pred_OLS,y_test_pred_causal_2SLS = OLSand2SLS(TrainVariants,TrainHeight,
                                                             TestVariants,TestHeight,
                                                             WholeTrainVariants,WholeTrainHeight)
            #Calculate the R2      
        R2ValuesOLS = getR2(TestHeight,y_test_pred_OLS,R2ValuesOLS)
            #Calculate the MAE      
        MAEValuesOLS = getMAE(TestHeight,y_test_pred_OLS,MAEValuesOLS)
            #Calculate the MSE      
        MSEValuesOLS = getMSE(TestHeight,y_test_pred_OLS,MSEValuesOLS)
            #Calculate the R2 Causal     
        R2Values2SLS = getR2(TestHeight,y_test_pred_causal_2SLS,R2Values2SLS)
            #Calculate the MAE Causal     
        MAEValues2SLS = getMAE(TestHeight,y_test_pred_causal_2SLS,MAEValues2SLS)
            #Calculate the MSE Causal     
        MSEValues2SLS = getMSE(TestHeight,y_test_pred_causal_2SLS,MSEValues2SLS)
            #Srtoring final errors 
        
    Results = {'Number of Features': numberoffeatures,
                     'OLS_R2': round(np.mean(R2ValuesOLS),7),'STD_OLS_R2': round(np.std(R2ValuesOLS),7),
                     '2SLS_R2': round(np.mean(R2Values2SLS),7),'STD_2SLS_R2': round(np.std(R2Values2SLS),7),
                     'RF_R2': round(np.mean(R2ValuesRF),7),'STD_RF_R2': round(np.std(R2ValuesRF),7),
                     '2SRF_R2': round(np.mean(R2Values2SRF),7),'STD_2SRF_R2': round(np.std(R2Values2SRF),7),
                     'KNN_R2': round(np.mean(R2ValuesKNN),7),'STD_KNN_R2': round(np.std(R2ValuesKNN),7),
                     '2SKNN_R2': round(np.mean(R2Values2SKNN),7),'STD_2SKNN_R2': round(np.std(R2Values2SKNN),7),
                     'OLS_MAE':  round(np.mean(MAEValuesOLS),7),'STD_OLS_MAE': round(np.std(MAEValuesOLS),7),
                     '2SLS_MAE': round(np.mean(MAEValues2SLS),7),'STD_2SLS_MAE': round(np.std(MAEValues2SLS),7),
                     'RF_MAE': round(np.mean(MAEValuesRF),7),'STD_RF_MAE': round(np.std(MAEValues2SRF),7),
                     '2SRF_MAE': round(np.mean(MAEValues2SRF),7),'STD_2SRF_MAE': round(np.std(MAEValues2SRF),7),
                     'KNN_MAE': round(np.mean(MAEValuesKNN),7),'STD_KNN_MAE': round(np.std(MAEValuesKNN),7),
                     '2SKNN_MAE': round(np.mean(MAEValues2SKNN),7),'STD_2SKNN_MAE': round(np.std(MAEValues2SKNN),7),
                     'OLS_MSE': round(np.mean(MSEValuesOLS),7),'STD_OLS_MSE': round(np.std(MSEValuesOLS),7),
                     '2SLS_MSE': round(np.mean(MSEValues2SLS),7),'STD_2SLS_MSE': round(np.std(MSEValues2SLS),7),
                     'RF_MSE': round(np.mean(MSEValuesRF),7),'STD_RF_MSE': round(np.std(MSEValuesRF),7),
                     '2SRF_MSE': round(np.mean(MSEValues2SRF),7),'STD_2SRF_MSE': round(np.std(MSEValues2SRF),7),
                     'KNN_MSE': round(np.mean(MSEValuesKNN),7),'STD_KNN_MSE': round(np.std(MSEValuesKNN),7),
                     '2SKNN_MSE': round(np.mean(MSEValues2SKNN),7),'STD_2SKNN_MSE': round(np.std(MSEValues2SKNN),7)}


    df = pd.concat([df, pd.DataFrame(data=Results, index=[0])], ignore_index=True)
    df.to_csv('D:/Users/ajarrah1/Desktop/HeightGenes/WholeCodeWithCausalInferenceResults/errors.csv')


    R2_names = ["OLS", "2SLS", "RF", 
               "2SRF", "KNN", "2SKNN"]
    R2_x_pos = np.arange(len(R2_names))
    R2_CTEs = [round(np.mean(R2ValuesOLS),7),
               round(np.mean(R2Values2SLS),7),
               round(np.mean(R2ValuesRF),7),
               round(np.mean(R2Values2SRF),7),
               round(np.mean(R2ValuesKNN),7),
               round(np.mean(R2Values2SKNN),7)]

    R2_error = [round(np.std(R2ValuesOLS),7),
               round(np.std(R2Values2SLS),7),
               round(np.std(R2ValuesRF),7),
               round(np.std(R2Values2SRF),7),
               round(np.std(R2ValuesKNN),7),
               round(np.std(R2Values2SKNN),7)]

    R2_fig, R2_ax = plt.subplots()
    R2_ax.bar(R2_x_pos,R2_CTEs, yerr=R2_error, align="center", alpha=0.5, ecolor='black', capsize=10)
    R2_ax.set_ylabel('R2')
    R2_ax.set_xticks(R2_x_pos)
    R2_ax.set_xticklabels(R2_names)
    R2_ax.set_title("R2 for OLS, 2SLS, RF, 2SRF, KNN, and 2SKNN for RFImp with "+str(numberoffeatures)+" features")
    R2_ax.yaxis.grid(True)

    plt.tight_layout()
    plt.savefig('D:/Users/ajarrah1/Desktop/HeightGenes/WholeCodeWithCausalInferenceResults/'+str(numberoffeatures)+" features R2.jpg")
    #plt.show

    MAE_names = ["OLS", "2SLS", "RF",
               "2SRF", "KNN", "2SKNN"]
    MAE_x_pos = np.arange(len(MAE_names))
    MAE_CTEs = [round(np.mean(MAEValuesOLS),7),
               round(np.mean(MAEValues2SLS),7),
               round(np.mean(MAEValuesRF),7),
               round(np.mean(MAEValues2SRF),7),
               round(np.mean(MAEValuesKNN),7),
               round(np.mean(MAEValues2SKNN),7)]

    MAE_error = [round(np.std(MAEValuesOLS),7),
               round(np.std(MAEValues2SLS),7),
               round(np.std(MAEValuesRF),7),
               round(np.std(MAEValues2SRF),7),
               round(np.std(MAEValuesKNN),7),
               round(np.std(MAEValues2SKNN),7)]

    MAE_fig, MAE_ax = plt.subplots()
    MAE_ax.bar(MAE_x_pos,MAE_CTEs, yerr=MAE_error, align="center", alpha=0.5, ecolor='black', capsize=10)
    MAE_ax.set_ylabel('MAE[inch]')
    MAE_ax.set_xticks(MAE_x_pos)
    MAE_ax.set_xticklabels(MAE_names)
    MAE_ax.set_title("MAE for OLS, 2SLS, RF, 2SRF, KNN, and 2SKNN for RFImp with "+str(numberoffeatures)+" features")
    MAE_ax.yaxis.grid(True)

    plt.tight_layout()
    plt.savefig('D:/Users/ajarrah1/Desktop/HeightGenes/WholeCodeWithCausalInferenceResults/'+str(numberoffeatures)+" features MAE.jpg")
    #plt.show

    MSE_names = ["OLS", "2SLS", "RF",
               "2SRF", "KNN", "2SKNN"]
    MSE_x_pos = np.arange(len(MSE_names))
    MSE_CTEs = [round(np.mean(MSEValuesOLS),7),
               round(np.mean(MSEValues2SLS),7),
               round(np.mean(MSEValuesRF),7),
               round(np.mean(MSEValues2SRF),7),
               round(np.mean(MSEValuesKNN),7),
               round(np.mean(MSEValues2SKNN),7)]

    MSE_error = [round(np.std(MSEValuesOLS),7),
               round(np.std(MSEValues2SLS),7),
               round(np.std(MSEValuesRF),7),
               round(np.std(MSEValues2SRF),7),
               round(np.std(MSEValuesKNN),7),
               round(np.std(MSEValues2SKNN),7)]

    MSE_fig, MSE_ax = plt.subplots()
    MSE_ax.bar(MSE_x_pos,MSE_CTEs, yerr=MSE_error, align="center", alpha=0.5, ecolor='black', capsize=10)
    MSE_ax.set_ylabel('MSE[inche^2]')
    MSE_ax.set_xticks(MSE_x_pos)
    MSE_ax.set_xticklabels(MSE_names)
    MSE_ax.set_title("MSE for OLS, 2SLS, RF, 2SRF, KNN, and 2SKNN for RFImp with "+str(numberoffeatures)+" features")
    MSE_ax.yaxis.grid(True)

    plt.tight_layout()
    plt.savefig('D:/Users/ajarrah1/Desktop/HeightGenes/WholeCodeWithCausalInferenceResults/'+str(numberoffeatures)+" features MSE.jpg")
    #plt.show
