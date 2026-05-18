def get_length_on(random_list):
    retlist = []
    count = 1
    i = 1
    while i < len(random_list):
        if random_list[i] == random_list[i-1] + 1:
            count += 1
        else:
            retlist.append(count)
            count = 1
        i += 1
    retlist.append(count)
    return retlist


def get_length_off(random_list):
    sorted_lst = sorted(random_list)
    retlist = []
    count = 0
    for i in range(1, len(sorted_lst)):
        if sorted_lst[i] - sorted_lst[i-1] > 1:
            count = sorted_lst[i] - sorted_lst[i-1] - 1
            retlist.append(count)
    if len(sorted_lst) == 1:
        retlist.append(0)
    return retlist

# def get_length_off(random_list, valeur=4000):
#     sorted_lst = sorted(random_list)
#     retlist = []
#     count = 0
#     for i in range(1, len(sorted_lst)):
#         if sorted_lst[i] - sorted_lst[i-1] > 1:
#             count = sorted_lst[i] - sorted_lst[i-1] - 1
#             retlist.append(count)
#     last_num = sorted_lst[-1]
#     if last_num < valeur:
#         distance_to_valeur = valeur - last_num
#         retlist.append(distance_to_valeur)
#     return retlist
