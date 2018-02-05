from __future__ import division

import numpy as np
import itertools

"""Change support and confidence here"""
minimum_support_percentage = 50.0
minimum_support = minimum_support_percentage / 100
print "Support is set to be ",
print minimum_support_percentage

minimum_confidence_percentage = 70.0
minimum_confidence = minimum_confidence_percentage / 100
print "Confidence is set to be ",
print minimum_confidence_percentage


total_number = 0
final_result = dict()
association_rules = []

def apriori():

    #Read from the modified file
    filename = "tempfile.txt"
    file = np.genfromtxt(filename, dtype=str, delimiter="\t")
    file = np.matrix(file)
    (rows, cols) = file.shape
    #print (rows,cols)

    #Required columns
    gene_cols = cols - 1

    global total_number
    #Generating item sets of length 1

    item_set = dict()
    possible_sets = []
    with open(filename, 'r')  as f:
        for line in f:
            line = line.split()
            #print len(line)
            for k in range(gene_cols):
                key = line[k]
                #print key
                if key in item_set:
                    item_set[key] = item_set[key] + 1
                else:
                    item_set[key] = 1

    for item in item_set.keys():
        val = item_set[item]
        support_item = val / rows
        # print support_item
        if (support_item < minimum_support):
            del (item_set[item])
        else:
            item_set[item] = support_item

    print "Length 1 frequent item sets: ",
    print item_set

    print "Number of length-1 frequent item sets: ",
    print len(item_set)

    total_number = total_number + len(item_set)

    final_result.update(item_set)

    #Generating frequent item sets of length >= 2

    n = len(item_set)
    #n=2
    start = 2

    for i in range(start,n+1):

        #Generating itemsets of length = 2:
        #item_set ----> frequent sets of length 1
        if i == 2:
            c = itertools.combinations(item_set, 2)
            frequent_sets = dict()
            for each_comb in c:
                frequency = 0
                with open(filename, 'r') as f:
                    for line in f:
                        if (set(each_comb) <= set(line.split())):
                            frequency = frequency + 1
                frequent_sets[each_comb] = frequency

            for item in frequent_sets.keys():
                val = frequent_sets[item]
                support_item = val / rows
                # print support_item
                if (support_item < minimum_support):
                    del (frequent_sets[item])
                else:
                    frequent_sets[item] = support_item
            print "Length 2 frequent item sets: ",
            print frequent_sets
            print "Number of length-2 frequent item sets: ",
            print len(frequent_sets)
            total_number = total_number + len(frequent_sets)
            final_result.update(frequent_sets)

        else:
            #Length of item sets > 2
            keys = frequent_sets.keys()
            possible_sets = []

            for j in keys:
                for k in keys:
                    comb = j + k
                    comb = set(comb)
                    comb = sorted(comb)
                    comb = tuple(comb)
                    # print type(comb)
                    if len(comb) == i:
                        possible_sets.append(comb)

            result_set = dict()
            for each_set in possible_sets:
                frequency = 0
                with open(filename, 'r') as f:
                    for line in f:
                        if (set(each_set) <= set(line.split())):
                            frequency = frequency + 1
                result_set[each_set] = frequency

            for item in result_set.keys():
                val = result_set[item]
                support_item = val / rows
                # print support_item
                if (support_item < minimum_support):
                    del (result_set[item])
                else:
                    result_set[item] = support_item

            #Exit the loop if no frequent set is generated
            if len(result_set) == 0:
                i = n + 1
            else:
                print "Number of length-",
                print i,
                print " frequent itemsets: ",
                print len(result_set)
                total_number = total_number + len(result_set)
                frequent_sets = result_set
                final_result.update(result_set)

    print "Total number of frequent item sets generated: ",
    print total_number
    print "Frequent item set: ",
    print final_result

    #Generating association rules for the frequent item sets
    for each_result in final_result:

        #Generate rules only if length of frequent item set is >= 2
        if( type(each_result) == tuple):
            complete_set_support = final_result[each_result]

            #Generate subsets of the given frequent set
            for i in range(1, len(each_result)):
                for j in itertools.combinations(each_result, i):

                    body = j #type:tuple
                    head = set(each_result) - set(j)

                    if len(body) == 1:
                        body = body[0]
                        if body in final_result.keys():
                            subset_support = final_result[body]
                            confidence = complete_set_support / subset_support

                            if confidence >= minimum_confidence:
                                head = tuple(head)
                                rule = []
                                rule.append(body)
                                rule.append(head)
                                association_rules.append(rule)

                    #subset has length >= 2
                    else:
                        for p in itertools.permutations(body):
                            if p in final_result:
                                body = p
                                subset_support = final_result[body]
                                confidence = complete_set_support / subset_support

                                if confidence >= minimum_confidence:
                                    head = tuple(head)
                                    rule = []
                                    rule.append(body)
                                    rule.append(head)
                                    association_rules.append(rule)

    print "Number of association rules generated: ",
    print len(association_rules)

    print "Association rules generated:"
    for rule in association_rules:
        print rule[0],
        print "====>",
        print rule[1]

def template1(str1,str2,feature_list):

    choice = str1
    number = -1

    #Checking if str2 is ANY/NONE or number
    try:
        if(int(str2)):
            number = int(str2)
    except ValueError:
        check = str2

    rule_count = 0
    resultant_rules = []

    if choice == "RULE":
        if number == -1:
            #ANY/NONE
            if check == "ANY":
                for gene in feature_list:
                    for rule in association_rules:
                        body = rule[0]
                        head = rule[1]

                        head_str = str(head)

                        if type(body) == str:
                            body_str = body
                            if gene == body:
                                result = "".join((body_str,"->",head_str))
                                if not result in resultant_rules:
                                    resultant_rules.append(result)
                                    print body,
                                    print "====>",
                                    print head

                        else:
                            body_str = str(body)
                            body_list = list(body)
                            if gene in body_list:
                                result = "".join((body_str, "->", head_str))
                                if not result in resultant_rules:
                                    resultant_rules.append(result)
                                    print body,
                                    print "====>",
                                    print head

                        head_list = list(head)
                        if gene in head_list:
                            result = "".join((body_str, "->", head_str))
                            if not result in resultant_rules:
                                resultant_rules.append(result)
                                print body,
                                print "====>",
                                print head


            if check == "NONE":
                feature_list_set = set(feature_list)
                for rule in association_rules:
                    body = rule[0]
                    head = rule[1]
                    if type(body) == str:
                        body_str = body
                        body = (body,)
                    else:
                        body_str = str(body)

                    body_set = set(body)
                    head_set = set(head)

                    head_str = str(head)

                    if feature_list_set.isdisjoint(body_set) and feature_list_set.isdisjoint(head_set):
                        result = "".join((body_str, "->", head_str))
                        if not result in resultant_rules:
                            resultant_rules.append(result)
                            print body,
                            print "====>",
                            print head

        else:
            #check is a number
            rules_added = []
            if (number <= len(feature_list)):
                # Generate combinations using number
                number = 1
                combinations = itertools.combinations(feature_list, number)
                for each_comb in combinations:
                    #print each_comb
                    #print type(each_comb) Ans: tuple
                    each_comb_set = set(each_comb)
                    #print each_comb_set

                    for rule in association_rules:
                        if not rule in rules_added:
                            body = rule[0]
                            head = rule[1]
                            head_set = set(head)
                            if type(body) == str:
                                body_set = (body,)
                                body_set = set(body_set)
                            else:
                                body_set = set(body)
                            if each_comb_set <= body_set:
                                rules_added.append(rule)

                            else:
                                if each_comb_set <= head_set:
                                    rules_added.append(rule)

                if len(feature_list) > 1:
                    # Checking for intersection
                    for rule in association_rules:
                        body = rule[0]
                        if type(body) == str:
                            body_tuple = (body,)
                        head = rule[1]
                        complete_rule = body_tuple + head
                        complete_rule = set(complete_rule)

                        feature_list_set = set(feature_list)

                        if feature_list_set <= complete_rule:
                            if rule in rules_added:
                                rules_added.remove(rule)

                #Adding rules to resultant_rules
                for rule in rules_added:
                    body = rule[0]
                    head = rule[1]
                    if type(body) == str:
                        body_str = body
                    else:
                        body_str = str(body)
                    head_str = str(head)
                    result = "".join((body_str, "->", head_str))
                    resultant_rules.append(result)
                    print body,
                    print "====>",
                    print head

    if choice == "BODY":

        if number == -1:
            #ANY/NONE
            if check == "ANY":
                for gene in feature_list:
                    for rule in association_rules:
                        body = rule[0]
                        head = rule[1]

                        head_str = str(head)

                        if type(body) == str:
                            body_str = body
                            if gene == body:
                                result = "".join((body_str,"->",head_str))
                                if not result in resultant_rules:
                                    resultant_rules.append(result)
                                    print body,
                                    print "====>",
                                    print head

                        else:
                            body_str = str(body)
                            body_list = list(body)
                            if gene in body_list:
                                result = "".join((body_str, "->", head_str))
                                if not result in resultant_rules:
                                    resultant_rules.append(result)
                                    print body,
                                    print "====>",
                                    print head

            if check == "NONE":
                feature_list_set = set(feature_list)
                for rule in association_rules:
                    body = rule[0]
                    head = rule[1]
                    if type(body) == str:
                        body_str = body
                        body = (body,)

                    else:
                        body_str = str(body)
                    body_set = set(body)

                    head_str = str(head)

                    if feature_list_set.isdisjoint(body_set):
                        result = "".join((body_str, "->", head_str))
                        if not result in resultant_rules:
                            resultant_rules.append(result)
                            print body,
                            print "====>",
                            print head


        else:
            #check is a number
            rules_added = []
            if (number <= len(feature_list)):
                # Generate combinations using number
                number = 1
                combinations = itertools.combinations(feature_list, number)
                for each_comb in combinations:
                    each_comb_set = set(each_comb)

                    for rule in association_rules:
                        if not rule in rules_added:
                            body = rule[0]
                            if type(body) == str:
                                body_set = (body,)
                                body_set = set(body_set)
                            else:
                                body_set = set(body)
                            if each_comb_set <= body_set:
                                rules_added.append(rule)

                if len(feature_list) > 1:
                    # Checking for intersection
                    for rule in association_rules:
                        body = rule[0]
                        if type(body) == str:
                            body_tuple = (body,)
                        else:
                            body_tuple = body
                        head = rule[1]
                        complete_rule = body_tuple
                        complete_rule = set(complete_rule)

                        feature_list_set = set(feature_list)

                        if feature_list_set <= complete_rule:
                            if rule in rules_added:
                                rules_added.remove(rule)

                # Adding rules to resultant_rules
                for rule in rules_added:
                    body = rule[0]
                    head = rule[1]
                    if type(body) == str:
                        body_str = body
                    else:
                        body_str = str(body)
                    head_str = str(head)
                    result = "".join((body_str, "->", head_str))
                    resultant_rules.append(result)
                    print body,
                    print "====>",
                    print head

    if choice == "HEAD":
        if number == -1:
            #ANY/NONE
            if check == "ANY":
                for gene in feature_list:
                    for rule in association_rules:
                        body = rule[0]
                        head = rule[1]

                        head_str = str(head)

                        if type(body) == str:
                            body_str = body
                        else:
                            body_str = str(body)

                        head_list = list(head)
                        if gene in head_list:
                            result = "".join((body_str, "->", head_str))
                            if not result in resultant_rules:
                                resultant_rules.append(result)
                                print body,
                                print "====>",
                                print head

            if check == "NONE":
                feature_list_set = set(feature_list)
                for rule in association_rules:
                    body = rule[0]
                    head = rule[1]
                    if type(body) == str:
                        body_str = body
                    else:
                        body_str = str(body)

                    head_set = set(head)

                    head_str = str(head)

                    if feature_list_set.isdisjoint(head_set):
                        result = "".join((body_str, "->", head_str))
                        if not result in resultant_rules:
                            resultant_rules.append(result)
                            print body,
                            print "====>",
                            print head

        else:
            # check is a number
            rules_added = []
            if (number <= len(feature_list)):
                # Generate combinations using number
                number = 1
                combinations = itertools.combinations(feature_list, number)
                for each_comb in combinations:
                    each_comb_set = set(each_comb)

                    for rule in association_rules:
                        if not rule in rules_added:
                            body = rule[0]
                            head = rule[1]
                            head_set = set(head)

                            if each_comb_set <= head_set:
                                rules_added.append(rule)

                if len(feature_list) > 1:
                    # Checking for intersection
                    for rule in association_rules:
                        head = rule[1]
                        complete_rule = head
                        complete_rule = set(complete_rule)

                        feature_list_set = set(feature_list)

                        if feature_list_set <= complete_rule:
                            if rule in rules_added:
                                rules_added.remove(rule)

                # Adding rules to resultant_rules
                for rule in rules_added:
                    body = rule[0]
                    head = rule[1]
                    if type(body) == str:
                        body_str = body
                    else:
                        body_str = str(body)
                    head_str = str(head)
                    result = "".join((body_str, "->", head_str))
                    resultant_rules.append(result)
                    print body,
                    print "====>",
                    print head

    print "No of rules: ",
    print len(resultant_rules)

    return resultant_rules


def template2(str1,number):

    choice = str1
    size = number

    resultant_rules =[]

    if choice == "RULE":
        for rule in association_rules:
            body = rule[0]
            head = rule[1]

            if type(body) == str:
                rule_size = 1 + len(list(head))
                body_str = body
            else:
                rule_size = len(list(body)) + len(list(head))
                body_str = str(body)
            head_str = str(head)

            if rule_size >= size:
                result = "".join((body_str, "->", head_str))
                if not result in resultant_rules:
                    resultant_rules.append(result)
                    print body,
                    print "====>",
                    print head


    if choice == "BODY":
        for rule in association_rules:
            body = rule[0]
            head = rule[1]
            if type(body) == str:
                body_size = 1
                body_str = body
            else:
                body_size = len(list(body))
                body_str = str(body)
            head_str = str(head)

            if body_size >= size:
                result = "".join((body_str, "->", head_str))
                if not result in resultant_rules:
                    resultant_rules.append(result)
                    print body,
                    print "====>",
                    print head


    if choice == "HEAD":
        for rule in association_rules:
            body = rule[0]
            head = rule[1]

            if type(body) == str:
                body_str = body
            else:
                body_str = str(body)

            head_str = str(head)
            head_size = len(list(head))

            if head_size >= size:
                result = "".join((body_str, "->", head_str))
                if not result in resultant_rules:
                    resultant_rules.append(result)
                    print body,
                    print "====>",
                    print head

    print "No of rules: ",
    print len(resultant_rules)

    return resultant_rules


def template3(*args):

    #Both the queries are from template 1
    if(len(args) == 7):
        print "Template 1"
        selection = args[0]
        choice1 = args[1]
        check1 = args[2]
        list1 = args[3]
        choice2 = args[4]
        check2 = args[5]
        list2 = args[6]

        result1 = template1(choice1,check1,list1)
        result2 = template1(choice2,check2,list2)

        result_set1 = set(result1)
        result_set2 = set(result2)

        if("or" in selection):
            result = result_set1.union(result_set2)
        else:
            result = result_set1.intersection(result_set2)

        print "Total number of rules: "
        print len(result)

        for i in result:
            print i


    #Query from both the templates
    elif(len(args) == 6):
        print "Template 1 and Template 2"
        selection = args[0]
        choice1 = args[1]
        check1 = args[2]
        list1 = args[3]

        choice2 = args[4]
        size2 = args[5]

        result1 = template1(choice1,check1,list1)
        result2 = template2(choice2,size2)

        result_set1 = set(result1)
        result_set2 = set(result2)

        if ("or" in selection):
            result = result_set1.union(result_set2)
        else:
            result = result_set1.intersection(result_set2)


        print "Total number of rules: "
        print len(result)

        for i in result:
            print i

    # Both the queries are from template 2
    elif(len(args) == 5):
        print "Template 2"
        selection = args[0]
        choice1 = args[1]
        size1 = args[2]
        choice2 = args[3]
        size2 = args[4]

        result1 = template2(choice1,size1)
        result2 = template2(choice2,size2)

        result_set1 = set(result1)
        result_set2 = set(result2)

        if ("or" in selection):
            result = result_set1.union(result_set2)
        else:
            result = result_set1.intersection(result_set2)


        print "Total number of rules: "
        print len(result)
        for i in result:
            print i


    else:
        print "Invalid query"



if __name__ == "__main__":
    apriori()
    print ""
    print "Executing Queries: "
    #template1("RULE", "ANY", ["G59_Up"])
    #template1("RULE", "NONE", ["G59_Up"])
    #template1("RULE", "1", ["G59_Up", "G10_Down"])

    #template1("BODY", "ANY", ["G59_Up"])
    #template1("BODY", "NONE", ["G59_Up"])
    #template1("BODY", "1", ["G59_Up", "G10_Down"])

    #template1("HEAD", "ANY", ["G59_Up"])
    #template1("HEAD", "NONE", ["G59_Up"])
    #template1("HEAD", "1", ["G59_Up", "G10_Down"])

    #template2("RULE", 3)
    #template2("BODY", 2)
    #template2("HEAD", 1)

    #template3("1or1", "BODY", "ANY", ["G10_Down"], "HEAD", "1", ["G59_Up"])
    #template3("1and1", "BODY", "ANY", ["G10_Down"], "HEAD", "1", ["G59_Up"])
    #template3("1or2", "BODY", "ANY", ["G10_Down"], "HEAD", 2)
    #template3("1and2", "BODY", "ANY", ["G10_Down"], "HEAD", 2)
    #template3("2or2", "BODY", 1, "HEAD", 2)
    #template3("2and2", "BODY", 1, "HEAD", 2)

    #template1("BODY","1",["G1_Up","G10_Down"])
    #template1("HEAD", "ANY", ["G1_Up"])
    #template2("HEAD",2)
    template3("2and2","BODY",1,"HEAD",2)
