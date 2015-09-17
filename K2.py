from __future__ import division;
import decimal;
import copy;
import math;
import random;
        
def readfile(filename):
    f = open(filename, "r");
    #read all lines;
    lines = f.readlines();
    wordstart = 0;
    word = "";
    nodedic = {};
    nodename = [];
    #read the first line into nodename[] as a list of each node's name;
    for char in lines[0]:
        # if the cursor hit a '"', then set wordstart to 1, prepare to record the word
        if char == '"' and wordstart == 0:
            wordstart = 1;
        elif char == '"' and wordstart == 1:        # if cursor hit a '"' again, then set the wordstart to 0, prepare to append the word between these two '"' to the list;
            wordstart = 0;
            nodedic[word[1:]] = 0;
            nodename.append(word[1:]);
            word = "";                              #set the word to "", prepare to record new word;
        if wordstart == 1 and char != "\x00" and char != "\r":          #if the wordstart is 1, means the cursor is inside two '"', then cancatenate the letter to word;
            word += char;
    database = [];
    #loop through all the other lines, build a tempdic for each line to represent the data for all columns in dataset
    #append the tempdic to the database, a list of all the training data;
    for line in lines[1:]:
        index = 0;
        tempdic = dict(nodedic);
        for char in line:
            if char != " " and char != "\n" and char != "\x00" and char != "\r":
                tempdic[nodename[index]] = int(char);
                index += 1;
        #del tempdic["INT"];
        database.append(tempdic);
    
    #nodename.remove("INT");
    return database, nodename;
#calculate the Nijk term in the equation;
def Nijk(i,parents, database): 
    result = [0,0,0];
    if parents == []:
        for entry in database:
            if entry[i] == 1:
                result[0] += 1;
            elif entry[i] == 2:
                result[1] += 1;
            elif entry [i] == 3:
                result[2] += 1;
        return result;
    else:
        newdatabase = copy.deepcopy(database);
        keylist = sorted(parents);
        #print "size: ",len(newdatabase);
        for entry in database:
            for key in keylist:
                if entry[key] != parents[key]:
                    newdatabase.remove(entry);
                    break;
        for entry in newdatabase:
            if entry[i] == 1:
                result[0] += 1;
            elif entry[i] == 2:
                result[1] += 1;
            elif entry [i] == 3:
                result[2] += 1;
        
        return result;
def multifactorial(res):        #calculate the last term in equation
    mfresult = 1;
    for k in res:
        mfresult *= math.factorial(k);
    return mfresult;
def sumN(res):              #calculate the sum of N
    sumresult = 0;
    for k in res:
        sumresult += k;
    return sumresult;

def gfunction(nodei,parenti,database):  # this is the g() function which is also the equation in the paper
    parentslist = [];
    parentdic = {};
    ri = 3;
    valist = [1,2,3];

    if len(parenti) == 0:               # if there's no parents, simply calculate the result once instead of multiply them together
        result = Nijk(nodei,[],database);
        mfresult = multifactorial(result);
        sumresult = sumN(result);
        gresult = mfresult*(math.factorial(ri-1)/decimal.Decimal(math.factorial(sumresult+ri-1)));
        #print "gresult: ", gresult;
        return gresult;
    
    #print "total list: ", valtotallist;
    while len(parentslist) < 3**len(parenti):       # if there are k parents, then each parent can have either of 3 values; so there are 3^k combinations
        for parentindex in range(0,len(parenti)):
            parentdic[parenti[parentindex]] = random.choice(valist);
        testnum = 1;                    #if this combination is already in parentdic then testnum = 0, try to get another new combination;
        for pl in parentslist:
            if len(parenti) == 2:
                if pl[parenti[0]] == parentdic[parenti[0]]  and pl[parenti[1]] == parentdic[parenti[1]]:
                    testnum = 0;    
                    break;
            elif len(parenti) == 1:
                if pl[parenti[0]] == parentdic[parenti[0]]:
                    testnum = 0;
                    break;
            elif len(parenti) == 3:
                if pl[parenti[0]] == parentdic[parenti[0]]  and pl[parenti[1]] == parentdic[parenti[1]] and pl[parenti[2]] == parentdic[parenti[2]]:
                    testnum = 0;
                    break;
        if testnum == 1:
            parentslist.append(dict(parentdic));
    gresult = 1;
    for parents in parentslist:
        result = Nijk(nodei,parents,database);
        mfresult = multifactorial(result);
        sumresult = sumN(result);
        #print "N:{} mf:{} sum:{}".format(result, mfresult, sumresult);
        #print "mfresult:{}, math.factorial(sumresult+ri-1):{}".format(mfresult, math.factorial(sumresult+ri-1));
        #print "mfresult/math.factorial(sumresult+ri-1):{}".format(float(mfresult/math.factorial(sumresult+ri-1)));
        gresult *= mfresult*(math.factorial(ri-1)/decimal.Decimal(math.factorial(sumresult+ri-1)));
    #print "len: {} result:{}".format(len(parentslist),gresult);
    return gresult;
    
# get the all the nodes before nodei in the initial nodeorder
def pred(nodei,nodeorder):
    predlist = []
    for nodename in nodeorder:
        if nodename == nodei:
            break;
        predlist.append(nodename);
    return predlist;
    
def K2(nodesets,nodeorder, u, database):
    
    for i in nodesets:
        pi = [];
        pold = gfunction(i,pi,database);
        OKToProceed = True;
        while OKToProceed and len(pi) < u:
            pnew = 0;
            predlist = pred(i, nodeorder);
            for p in pi:
                predlist.remove(p);
            for z in predlist:
                newpi = list(pi);
                newpi.append(z)
                #print "newpi: ", newpi;
                pnewprime = gfunction(i, newpi,database);
                if pnewprime > pnew:
                    pnew = pnewprime;
                    znew = z;
            if pnew > pold:
                pold = pnew;
                pi.append(znew);
            else:
                OKToProceed = False;
        print "Node: {}, Parents of this node: {}".format(i, pi);
            
    
    
#this is the main function where this script begins to execute
if __name__ == "__main__":
    filename = "sachs_intervention_honghzu.txt";
    database, nodename = readfile(filename);
    #nodeorder = ["PKC","plcg","praf","PKA","PIP3","PIP2","pmek","p44.42","pjnk","pakts473","P38"];
    nodeorder = ["PKC","plcg","PKA","PIP3","P38","praf","pjnk","PIP2","pmek","p44.42","pakts473"];
    #nodeorder = ["PKC","PKA","P38","praf","pjnk","pmek","p44.42","pakts473","plcg","PIP3","PIP2"];
    K2(nodename, nodeorder, 3, database);
    print len(database);
    
