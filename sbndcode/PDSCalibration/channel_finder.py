# Code that takes the argument of PMT channel number and gives you corresponding index to use with SPE finding code
# H. Parkinson for SBND, 2024


import sys

pmt_numbers = [17, 71, 41, 15, 69, 13, 67, 65, 16, 70, 40, 14, 68, 12, 66, 39, 11, 9, 63, 37, 7, 61, 64, 38, 10, 8, 62, 36, 6, 60, 95, 149, 119, 93, 147, 91, 145, 143, 94, 148, 118, 92, 146, 90, 144, 117, 89, 87, 141, 115, 85, 139, 142, 116, 88, 86, 140, 114, 84, 138, 173, 227, 197, 171, 225, 169, 223, 221, 172, 226, 196, 170, 224, 168, 222, 195, 167, 165, 219, 193, 163, 217, 220, 194, 166, 164, 218, 192, 162, 216, 251, 305, 275, 249, 303, 247, 301, 299, 250, 304, 274, 248, 302, 246, 300, 273, 245, 243, 297, 271, 241, 295, 298, 272, 244, 242, 296, 270, 240, 294]


sorted_pmt_numbers = sorted(pmt_numbers)

def find_index(PMT):
	print(sorted_pmt_numbers.index(PMT))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python greeting.py <PMT>")
        sys.exit(1)
    
    PMT = int(sys.argv[1])
    
    if PMT not in sorted_pmt_numbers:
        print("That number isn't on the list of PMT channels!")
        sys.exit(1)



    find_index(PMT)













