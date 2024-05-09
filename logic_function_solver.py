# coding: utf-8

"""
@Author: JC
"""

from texttable import Texttable


def minterm_to_binary(minterm, num_vars):
    """
    Convert a minterm to its binary representation.
    """
    binary = bin(minterm)[2:].zfill(num_vars)
    return binary


def get_ones_count(binary):
    """
    Return the number of '1' in a binary term.
    """
    return binary.count('1')


def get_dashes_count(binary):
    """
    Return the number of '-' in a binary term.
    """
    return binary.count('-')


def combine_terms(terms, one_counter, minterms_set):
    """
    Combine and group terms with one difference
    """
    combined_terms = []
    used_terms = []
    terms_set = []
    used_terms_set = []

    for i in range(len(terms)):
        for j in range(i + 1, len(terms)):
            term1 = terms[i]
            term2 = terms[j]
            term1_set = minterms_set[i]
            term2_set = minterms_set[j]
            if get_ones_count(term1) != get_ones_count(term2) - 1:  # Combines terms with only one '1' difference
                continue
            temp = 0
            for k in range(len(term1)):
                if term1[k] == "-" and term2[k] == "-":  # only terms with same dash position and quantity will be considered
                    temp += 1
            if temp != one_counter:
                continue

            combined_term = ''
            for k in range(len(term1)):
                if term1[k] == term2[k]:
                    combined_term += term1[k]
                else:
                    combined_term += '-'
            if get_dashes_count(combined_term) == (one_counter + 1):  # the number of '-' must be equal to the generation of iterations
                combined_terms.append(combined_term)
                used_terms.append(term1)
                used_terms.append(term2)
                used_terms_set.append(tuple(term1_set))
                used_terms_set.append(tuple(term2_set))
                temp = list(tuple(term1_set) + tuple(term2_set))
                temp.sort()
                terms_set.append(tuple(temp))

    if used_terms:
        for term in terms:
            if term not in used_terms:
                combined_terms.append(term)  # the term can't be combined with all terms will become PI

    if used_terms_set:
        for term in minterms_set:
            if term not in used_terms_set:
                terms_set.append(term)

    if not combined_terms:  # end of iteration and return a list of PI
        return True, terms, minterms_set
    return False, combined_terms, terms_set


def check_only(minterm, terms_set, current_set):
    """
    Determine if a minterm is in other minterm groups.
    """
    for term_set in terms_set:
        if term_set == current_set:
            continue
        if minterm in term_set:
            return False
    return True


def remove_covered_by_epi(pi_terms_set, epi_terms_set, pi, dc):
    """
    If a minterm is covered by EPI, remove it from its own PI group
    """
    pi_terms = [i for i in pi_terms_set if i not in epi_terms_set]
    epi_temp = []
    epi_set_temp = []
    perfect_set = []
    for term_set in pi_terms:
        set_temp = ()
        for j in term_set:
            if check_only(j, epi_terms_set, term_set):
                if j not in dc:
                    set_temp += (j,)
        if not set_temp == ():
            perfect_set.append(term_set)
            epi_temp.append(pi[pi_terms_set.index(term_set)])
            epi_set_temp.append(set_temp)

    return epi_temp, epi_set_temp, perfect_set


def find_prime_implicants(num_vars, minterms, minterms_set):
    """
    Find the prime implicants using the tabulation method.
    """
    terms = [minterm_to_binary(minterm, num_vars) for minterm in minterms]
    # print(terms)
    for i in range(num_vars):
        done, terms, minterms_set = combine_terms(terms, i, minterms_set)
        # print(terms)
        # print(minterms_set)
        if done:
            break
    return terms, minterms_set


def find_epi(terms_set, pi, dc):
    """
    Find EPI
    """
    epi_terms_set = []
    epi_terms = []
    for term_set in terms_set:
        for j in term_set:
            if check_only(j, terms_set, term_set):
                if j not in dc:  # exclude the don't care
                    epi_terms_set.append(term_set)
                    epi_terms.append(pi[terms_set.index(term_set)])
                break

    pi_terms_set = [i for i in terms_set if i not in epi_terms_set]
    pi_terms = [i for i in pi if i not in epi_terms]

    while True:
        perfect_term, removed_epi_set, perfect_set = remove_covered_by_epi(pi_terms_set, epi_terms_set, pi_terms, dc)
        max_len = -1
        epi_term = ""
        epi_set = ()
        # most covered PI is an EPI
        for term in removed_epi_set:
            if len(term) > max_len:
                max_len = len(term)
                epi_set = perfect_set[removed_epi_set.index(term)]
                epi_term = perfect_term[removed_epi_set.index(term)]

        if not epi_term:
            break
        else:
            epi_terms.append(epi_term)
            epi_terms_set.append(epi_set)

    #  remove redundant term
    for term_set in epi_terms_set:
        count = 0
        for j in term_set:
            if not check_only(j, epi_terms_set, term_set):
                count += 1
            else:
                if j in dc:
                    count += 1
        if count == len(term_set):
            epi_terms.pop(epi_terms_set.index(term_set))
            epi_terms_set.remove(term_set)

    return epi_terms, epi_terms_set


def print_table(num_vars, prime_implicants, pi_set, epi_terms, epi_set, MINTERM):
    is_epi = ["^_^" if prime_implicants[i] in epi_terms else "" for i in range(len(prime_implicants))]
    symbol_groups = []
    plus = False
    for term in prime_implicants:
        symbol = ""
        for j in range(num_vars):
            if term[j] != '-':
                if not MINTERM:
                    if plus:
                        symbol += " + "
                    else:
                        plus = True
                symbol += chr(ord('A') + j)

            if term[j] == '0':
                if MINTERM:
                    symbol += "'"
            if term[j] == '1':
                if not MINTERM:
                    symbol += "'"
        plus = False

        symbol_groups.append(symbol)

    table = Texttable()
    table.set_cols_align(["c", "c", "c", "c"])
    table.set_header_align(["c", "c", "c", "c"])
    header = ["Prime Implicant", "Prime Implicant Group", "Is EPI", "Symbol(s)"]
    table.add_row(header)
    combined = list([prime_implicants[i], pi_set[i], is_epi[i], symbol_groups[i]] for i in range(len(prime_implicants)))
    for i in combined:
        table.add_row(i)

    print(table.draw())
    print()
    symbol_groups = []
    for term in epi_terms:
        symbol = ""
        for j in range(num_vars):
            if term[j] != '-':
                if not MINTERM:
                    if plus:
                        symbol += " + "
                    else:
                        plus = True
                symbol += chr(ord('A') + j)

            if term[j] == '0':
                if MINTERM:
                    symbol += "'"
            if term[j] == '1':
                if not MINTERM:
                    symbol += "'"
        plus = False

        symbol_groups.append(symbol)

    table = Texttable()
    table.set_cols_align(["c", "c", "c"])
    table.set_header_align(["c", "c", "c"])
    table.header(["Essential Prime Implicant", "Essential Prime Implicant Group", "Symbol(s)"])
    combined = list([epi_terms[i], epi_set[i], symbol_groups[i]] for i in range(len(epi_terms)))
    for i in combined:
        table.add_row(i)
    print(table.draw())

    print("F(", end="")
    for i in range(num_vars):
        print(chr(ord('A') + i), end="")
        if i != num_vars - 1:
            print(",", end="")

    print("): ", end="")

    plus = False
    for i in symbol_groups:
        if MINTERM:
            if plus:
                print(" + ", end="")
            else:
                plus = True
            print(i, end="")
        else:
            print(f"({i})", end="")

    print()


def main():
    num_vars = int(input("Enter the number of variables: "))
    minterms = list(map(int, input("Enter the minterms separated by spaces: ").split()))
    dc = list(map(int, input("Enter the don't cares(otherwise leave it blank): ").split()))

    """
    For testing purpose
    """
    # num_vars = 5
    # minterms = list(map(int, "0 4 8 9 10 11 12 13 15".split()))
    # minterms = list(map(int, "0 1 2 3 4 6 7 11 12 15".split()))
    # minterms = list(map(int, "0 1 2 4 5 6 8 9 10 14 15".split()))
    # minterms = list(map(int, "0 1 2 4 5 9 11 12 15".split()))
    # minterms = list(map(int, "0 1 2 4 5 9 10 11 12 15".split()))
    # minterms = list(map(int, "0 1 2 4 5 7 8 9 11 12 15 19 20 22 26 27 28 30 34".split()))
    # minterms = list(map(int, "1 3 6 8 9 29 30 31 33 37 38 40 41 43 44 47 48 51 55 57 59 60 61 62 63".split()))
    # dc = [3, 14]

    minterms_set = []
    maxterms = []
    maxterms_set = []

    if dc:
        for i in dc:
            minterms.append(i)
            maxterms.append(i)

    minterms.sort()
    minterms_set.extend(list((i,) for i in minterms))
    maxterms.extend(list(i for i in range(2 ** num_vars) if i not in minterms))  # find maxterms from given minterms
    maxterms.sort()
    maxterms_set.extend(list((i,) for i in maxterms))

    prime_implicants, pi_set = find_prime_implicants(num_vars, minterms, minterms_set)
    prime_implicants = list(dict.fromkeys(prime_implicants))  # remove minterms duplicates
    pi_set = list(dict.fromkeys(pi_set))  # remove minterms set duplicates
    epi_terms, epi_set = find_epi(pi_set, prime_implicants, dc)
    epi_terms = list(dict.fromkeys(epi_terms))  # remove epi terms duplicates
    epi_set = list(dict.fromkeys(epi_set))  # remove epi set terms duplicates

    print("\nThe SOP expression of the logic function:")
    print_table(num_vars, prime_implicants, pi_set, epi_terms, epi_set, True)
    print()
    print("_"*76)
    print()
    print("The POS expression of the logic function:")
    prime_implicants, pi_set = find_prime_implicants(num_vars, maxterms, maxterms_set)
    prime_implicants = list(dict.fromkeys(prime_implicants))  # remove maxterms duplicates
    pi_set = list(dict.fromkeys(pi_set))  # remove maxterms_set duplicates
    epi_terms, epi_set = find_epi(pi_set, prime_implicants, dc)
    epi_terms = list(dict.fromkeys(epi_terms))  # remove epi terms duplicates
    epi_set = list(dict.fromkeys(epi_set))  # remove epi set terms duplicates

    print_table(num_vars, prime_implicants, pi_set, epi_terms, epi_set, False)


if __name__ == "__main__":
    main()
