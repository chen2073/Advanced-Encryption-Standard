from BitVector import *
from time import time
from x931 import x931
from copy import deepcopy


AES_modulus = BitVector(bitstring='100011011')
subBytesTable = [99, 124, 119, 123, 242, 107, 111, 197, 48, 1, 103, 43, 254, 215, 171, 118, 202, 130, 201, 125, 250, 89, 71, 240, 173, 212, 162, 175, 156, 164, 114, 192, 183, 253, 147, 38, 54, 63, 247, 204, 52, 165, 229, 241, 113, 216, 49, 21, 4, 199, 35, 195, 24, 150, 5, 154, 7, 18, 128, 226, 235, 39, 178, 117, 9, 131, 44, 26, 27, 110, 90, 160, 82, 59, 214, 179, 41, 227, 47, 132, 83, 209, 0, 237, 32, 252, 177, 91, 106, 203, 190, 57, 74, 76, 88, 207, 208, 239, 170, 251, 67, 77, 51, 133, 69, 249, 2, 127, 80, 60, 159, 168, 81, 163, 64, 143, 146, 157, 56, 245, 188, 182, 218, 33, 16, 255, 243, 210, 205, 12, 19, 236, 95, 151, 68, 23, 196, 167, 126, 61, 100, 93, 25, 115, 96, 129, 79, 220, 34, 42, 144, 136, 70, 238, 184, 20, 222, 94, 11, 219, 224, 50, 58, 10, 73, 6, 36, 92, 194, 211, 172, 98, 145, 149, 228, 121, 231, 200, 55, 109, 141, 213, 78, 169, 108, 86, 244, 234, 101, 122, 174, 8, 186, 120, 37, 46, 28, 166, 180, 198, 232, 221, 116, 31, 75, 189, 139, 138, 112, 62, 181, 102, 72, 3, 246, 14, 97, 53, 87, 185, 134, 193, 29, 158, 225, 248, 152, 17, 105, 217, 142, 148, 155, 30, 135, 233, 206, 85, 40, 223, 140, 161, 137, 13, 191, 230, 66, 104, 65, 153, 45, 15, 176, 84, 187, 22]

def gee(keyword, round_constant, byte_sub_table):
    # g function
    rotated_word = keyword.deep_copy()
    rotated_word << 8
    newword = BitVector(size = 0)
    for i in range(4):
        newword += BitVector(intVal = byte_sub_table[rotated_word[8*i:8*i+8].intValue()], size=8)
    newword[:8] ^= round_constant
    round_constant = round_constant.gf_multiply_modular(BitVector(intVal=0x02), AES_modulus, 8)
    return newword, round_constant


def gen_key_schedule_256(key_bv):
    global subBytesTable
    byte_sub_table = subBytesTable
    #  We need 60 keywords (each keyword consists of 32 bits) in the key schedule for
    #  256 bit AES. The 256-bit AES uses the first four keywords to xor the input
    #  block with.  Subsequently, each of the 14 rounds uses 4 keywords from the key
    #  schedule. We will store all 60 keywords in the following list:
    key_words = [None for i in range(60)]
    round_constant = BitVector(intVal = 0x01, size=8)
    for i in range(8):
        key_words[i] = key_bv[i*32 : i*32 + 32]
    for i in range(8,60):
        if i%8 == 0:
            kwd, round_constant = gee(key_words[i-1], round_constant, byte_sub_table)
            key_words[i] = key_words[i-8] ^ kwd
        elif (i - (i//8)*8) < 4:
            key_words[i] = key_words[i-8] ^ key_words[i-1]
        elif (i - (i//8)*8) == 4:
            key_words[i] = BitVector(size = 0)
            for j in range(4):
                key_words[i] += BitVector(intVal=byte_sub_table[key_words[i-1][8*j:8*j+8].intValue()], size = 8)
            key_words[i] ^= key_words[i-8]
        elif ((i - (i//8)*8) > 4) and ((i - (i//8)*8) < 8):
            key_words[i] = key_words[i-8] ^ key_words[i-1]
        else:
            sys.exit("error in key scheduling algo for i = %d" % i)
    return key_words


key_s = gen_key_schedule_256(BitVector(textstring="applesbananaspeachesstrawberries"))

def AES(bv, key):
    if bv.length() != 128:
        raise ValueError("block is not 128 bits")

    if not isinstance(key, BitVector):
        raise ValueError("key must be a BitVector")

    key_words = deepcopy(key_s)
    global subBytesTable
    statearray = [[0 for _ in range(4)] for _ in range(4)]
    temp_shift = [0] * 4

    # Initialize state array and add 1st round key
    for i in range(4):
        for j in range(4):
            statearray[i][j] = bv[32 * i + 8 * j:32 * i + 8 * (j + 1)]
            statearray[i][j] ^= key_words[i][8 * j:8 + (8 * j)]

    for rounds in range(14):
        # SubBytes
        for i in range(4):
            for j in range(4):
                statearray[i][j] = BitVector(intVal=subBytesTable[int(statearray[i][j])])

        # ShiftRows
        for i in range(1, 4):
            for j in range(0, 4):
                temp_shift[(j - i) % 4] = statearray[j][i]
            for j in range(0, 4):
                statearray[j][i] = temp_shift[j]

        # ColumnMixing
        if rounds != 13:
            two_times = BitVector(bitstring='00000010')
            three_times = BitVector(bitstring='00000011')
            for i in range(4):
                temp = (two_times.gf_multiply_modular(statearray[i][0], AES_modulus, 8)) ^ (
                    three_times.gf_multiply_modular(statearray[i][1], AES_modulus, 8)) ^ statearray[i][2] ^ \
                       statearray[i][3]
                temp1 = (two_times.gf_multiply_modular(statearray[i][1], AES_modulus, 8)) ^ (
                    three_times.gf_multiply_modular(statearray[i][2], AES_modulus, 8)) ^ statearray[i][3] ^ \
                        statearray[i][0]
                temp2 = (two_times.gf_multiply_modular(statearray[i][2], AES_modulus, 8)) ^ (
                    three_times.gf_multiply_modular(statearray[i][3], AES_modulus, 8)) ^ statearray[i][0] ^ \
                        statearray[i][1]
                temp3 = (two_times.gf_multiply_modular(statearray[i][3], AES_modulus, 8)) ^ (
                    three_times.gf_multiply_modular(statearray[i][0], AES_modulus, 8)) ^ statearray[i][1] ^ \
                        statearray[i][2]

                statearray[i][0] = temp
                statearray[i][1] = temp1
                statearray[i][2] = temp2
                statearray[i][3] = temp3

        # Add round key
        for i in range(4):
            for j in range(4):
                statearray[i][j] ^= key_words[(4 * (rounds + 1)) + i][8 * j:8 + (8 * j)]

        # transform 4x4 array to linear list
        output = ""
        for i in range(4):
            for j in range(4):
                # output is a binary file
                output += statearray[i][j].get_bitvector_in_hex()

    return BitVector(hexstring=output)


def ctr_aes_image(iv, image_file, out_file, key):
    # header = b""
    with open(image_file, "rb") as file, open("temp.ppm", "wb") as temp:
        lines = file.readlines()
        header = b"".join(lines[0:3])
        image = b"".join(lines[3:])
        temp.write(image)

    with open(out_file, "ab") as fout:
        fout.write(header)
        bv = BitVector(filename="temp.ppm")
        # initialize sufficient random number using x931
        dt = BitVector(intVal=int((time() * 1000000)), size=64) + BitVector(intVal=int((time() * 1000000)), size=64)
        nonce = x931(iv, dt, 1, key)[0]
        # print("nonce:", int(nonce))
        i=0
        while bv.more_to_read:
            bitvec = bv.read_bits_from_file(128)
            if len(bitvec) > 0:
                if len(bitvec) != 128:
                    bitvec.pad_from_right(128 - len(bitvec))
            try:
                cipher_block = AES(nonce, key) ^ bitvec
            except Exception as e:
                print("nonce length:", nonce.length())
                print("bitvec length:", bitvec.length())
                raise e
            cipher_block.write_to_file(fout)
            nonce = BitVector(intVal=(int(nonce) + 1), size=128)
            print(i)
            i += 1


if __name__ == "__main__":
    import os
    os.remove('enc_image.ppm') if os.path.isfile('enc_image.ppm') else print()
    v0 = BitVector(textstring="computersecurity")
    key = BitVector(textstring="applesbananaspeachesstrawberries")
    ctr_aes_image(v0, 'image.ppm', 'enc_image.ppm', key)
    print("done")

