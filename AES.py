#!/usr/bin/env python3

from BitVector import *
from functools import reduce
from copy import deepcopy

# global variables
AES_modulus = BitVector(bitstring='100011011')
subBytesTable = []  # for encryption
invSubBytesTable = []  # for decryption


def gen_subbytes_table():
    subBytesTable = []
    c = BitVector(bitstring='01100011')
    for i in range(0, 256):
        a = BitVector(intVal=i, size=8).gf_MI(AES_modulus, 8) if i != 0 else BitVector(intVal=0)
        a1,a2,a3,a4 = [a.deep_copy() for x in range(4)]
        a ^= (a1 >> 4) ^ (a2 >> 5) ^ (a3 >> 6) ^ (a4 >> 7) ^ c
        subBytesTable.append(int(a))
    return subBytesTable


def genTables():
    c = BitVector(bitstring='01100011')
    d = BitVector(bitstring='00000101')
    for i in range(0, 256):
        # For the encryption SBox
        a = BitVector(intVal = i, size=8).gf_MI(AES_modulus, 8) if i != 0 else BitVector(intVal=0)
        # For bit scrambling for the encryption SBox entries:
        a1,a2,a3,a4 = [a.deep_copy() for x in range(4)]
        a ^= (a1 >> 4) ^ (a2 >> 5) ^ (a3 >> 6) ^ (a4 >> 7) ^ c
        subBytesTable.append(int(a))
        # For the decryption Sbox:
        b = BitVector(intVal = i, size=8)
        # For bit scrambling for the decryption SBox entries:
        b1,b2,b3 = [b.deep_copy() for x in range(3)]
        b = (b1 >> 2) ^ (b2 >> 5) ^ (b3 >> 7) ^ d
        check = b.gf_MI(AES_modulus, 8)
        b = check if isinstance(check, BitVector) else 0
        invSubBytesTable.append(int(b))


def gen_round_keys(key_bv):
    key_words = gen_key_schedule_256(key_bv)
    num_rounds = 14
    round_keys = [None for i in range(num_rounds + 1)]
    for i in range(num_rounds + 1):
        round_keys[i] = (key_words[i * 4] + key_words[i * 4 + 1] + key_words[i * 4 + 2] + key_words[i * 4 + 3])
    return round_keys


def gen_key_schedule_256(key_bv):
    byte_sub_table = gen_subbytes_table()
    #  We need 60 keywords (each keyword consists of 32 bits) in the key schedule for
    #  256 bit AES. The 256-bit AES uses the first four keywords to xor the input
    #  block with.  Subsequently, each of the 14 rounds uses 4 keywords from the key
    #  schedule. We will store all 60 keywords in the following list:
    key_words = [None for i in range(60)]
    round_constant = BitVector(intVal=0x01, size=8)
    for i in range(8):
        key_words[i] = key_bv[i*32: i*32 + 32]
    for i in range(8,60):
        if i%8 == 0:
            kwd, round_constant = gee(key_words[i-1], round_constant, byte_sub_table)
            key_words[i] = key_words[i-8] ^ kwd
        elif (i - (i//8)*8) < 4:
            key_words[i] = key_words[i-8] ^ key_words[i-1]
        elif (i - (i//8)*8) == 4:
            key_words[i] = BitVector(size=0)
            for j in range(4):
                key_words[i] += BitVector(intVal=byte_sub_table[key_words[i-1][8*j:8*j+8].intValue()], size=8)
            key_words[i] ^= key_words[i-8]
        elif ((i - (i//8)*8) > 4) and ((i - (i//8)*8) < 8):
            key_words[i] = key_words[i-8] ^ key_words[i-1]
        else:
            sys.exit("error in key scheduling algorithm for i = %d" % i)
    return key_words


def gee(keyword, round_constant, byte_sub_table):
    # This is the g() function you see in Figure 4 of Lecture 8.
    rotated_word = keyword.deep_copy()
    rotated_word << 8
    newword = BitVector(size=0)
    for i in range(4):
        newword += BitVector(intVal=byte_sub_table[rotated_word[8*i:8*i+8].intValue()], size=8)
    newword[:8] ^= round_constant
    round_constant = round_constant.gf_multiply_modular(BitVector(intVal=0x02), AES_modulus, 8)
    return newword, round_constant


def SubByte(state):
    # return [[BitVector(size=8, intVal=subBytesTable[state[i][j].int_val()]) for j in range(4)] for i in range(4)]
    for i in range(4):
        for j in range(4):
            state[j][i] = BitVector(size=8, intVal=subBytesTable[state[i][j].int_val()])
    return state


def ShiftRows(state):
    state[1] = state[1][1:] + state[1][:1]
    state[2] = state[2][2:] + state[2][:2]
    state[3] = state[3][3:] + state[3][:3]
    return state


def MixColumns(state):
    temp = deepcopy(state)
    hex_02 = BitVector(hexstring='02')
    hex_03 = BitVector(hexstring='03')
    for i in range(4):
        state[0][i] = temp[0][i].gf_multiply_modular(hex_02, AES_modulus, 8) ^ temp[1][i].gf_multiply_modular(hex_03, AES_modulus, 8) ^ temp[2][i] ^ temp[3][i]
        state[1][i] = temp[0][i] ^ temp[1][i].gf_multiply_modular(hex_02, AES_modulus, 8) ^ temp[2][i].gf_multiply_modular(hex_03, AES_modulus, 8) ^ temp[3][i]
        state[2][i] = temp[0][i] ^ temp[1][i] ^ temp[2][i].gf_multiply_modular(hex_02, AES_modulus, 8) ^ temp[3][i].gf_multiply_modular(hex_03, AES_modulus, 8)
        state[3][i] = temp[0][i].gf_multiply_modular(hex_03, AES_modulus, 8) ^ temp[1][i] ^ temp[2][i] ^ temp[3][i].gf_multiply_modular(hex_02, AES_modulus, 8)

    return state


# 1st round adding round key
def Add1stRoundKey(block, word):
    return block ^ word


# state array and round key
def AddRoundKey(state, word):
    #print(state)
    word_matrix = [[word[32*i + 8*j: 32*i+8*(j+1)] for i in range(4)] for j in range(4)]
    #print(word_matrix)
    for i in range(4):
        for j in range(4):
            state[i][j] ^= word_matrix[i][j]
    return state


def encrypt():
    # populate s-table
    genTables()

    # key_bv = BitVector(filename="key.txt")
    key_bv = BitVector(textstring="anexaminedlifeistrulyworthliving")
    msg_bv = BitVector(filename="message.txt")

    # generate round keys: key size must be 256
    roundKeys = gen_round_keys(key_bv)

    output = BitVector(size=0)
    while msg_bv.more_to_read:
        # fill one block with 128 bits
        block = msg_bv.read_bits_from_file(128)
        # check padding needed
        if block.length() < 128:
            block.pad_from_right(128-block.length())

        # adding 1st round key
        block = Add1stRoundKey(block, roundKeys[0])
        # construct to 4x4 state array
        # state_array = [[block[32*i + 8*j: 32*i+8*(j+1)] for i in range(4)] for j in range(4)]
        state_array = [[0 for _ in range(4)] for _ in range(4)]
        for i in range(4):
            for j in range(4):
                state_array[j][i] = BitVector(size=8, intVal=block[32 * i + 8 * j:32 * i + 8 * (j + 1)].int_val())

        # 14 rounds of processing
        for i in range(1, 15):
            # last round
            if i == 14:
                state_array = ShiftRows(SubByte(state_array))
            else:
                state_array = MixColumns(ShiftRows(SubByte(state_array)))

            state_array = AddRoundKey(state_array, roundKeys[i])

        # convert 4x4 state array to one dimensional array
        for i in range(4):
            for j in range(4):
                output += state_array[i][j]

    # with open("encrypted.txt", "w") as file:
    #     file.write(output.get_bitvector_in_hex())
    fp = open("encrypted_b.txt", "wb")
    output.write_to_file(fp)
    fp.close()


def InvSubBytes(state):
    for i in range(4):
        for j in range(4):
                state[j][i] = BitVector(size=8, intVal=invSubBytesTable[state[i][j].int_val()])
    return state


def InvShiftRows(state):
    state[1] = state[1][-1:]+state[1][:-1]
    state[2] = state[2][-2:]+state[2][:-2]
    state[3] = state[3][-3:]+state[3][:-3]
    return state


def InvMixColumns(state):
    hex_0E = BitVector(hexstring='0e')
    hex_0B = BitVector(hexstring='0b')
    hex_0D = BitVector(hexstring='0d')
    hex_09 = BitVector(hexstring='09')
    temp = deepcopy(state)
    for i in range(4):
        state[0][i] = temp[0][i].gf_multiply_modular(hex_0E, AES_modulus, 8) ^ temp[1][i].gf_multiply_modular(hex_0B, AES_modulus, 8) ^ temp[2][i].gf_multiply_modular(hex_0D, AES_modulus, 8) ^ temp[3][i].gf_multiply_modular(hex_09, AES_modulus, 8)
        state[1][i] = temp[0][i].gf_multiply_modular(hex_09, AES_modulus, 8) ^ temp[1][i].gf_multiply_modular(hex_0E, AES_modulus, 8) ^ temp[2][i].gf_multiply_modular(hex_0B, AES_modulus, 8) ^ temp[3][i].gf_multiply_modular(hex_0D, AES_modulus, 8)
        state[2][i] = temp[0][i].gf_multiply_modular(hex_0D, AES_modulus, 8) ^ temp[1][i].gf_multiply_modular(hex_09, AES_modulus, 8) ^ temp[2][i].gf_multiply_modular(hex_0E, AES_modulus, 8) ^ temp[3][i].gf_multiply_modular(hex_0B, AES_modulus, 8)
        state[3][i] = temp[0][i].gf_multiply_modular(hex_0B, AES_modulus, 8) ^ temp[1][i].gf_multiply_modular(hex_0D, AES_modulus, 8) ^ temp[2][i].gf_multiply_modular(hex_09, AES_modulus, 8) ^ temp[3][i].gf_multiply_modular(hex_0E, AES_modulus, 8)

    return state


def decrypt():
    if subBytesTable == [] or invSubBytesTable == []:
        genTables()

    key_bv = BitVector(textstring="anexaminedlifeistrulyworthliving")
    # text = open("encrypted.txt", "r").read()
    # msg_bv = BitVector(hexstring=text)
    msg_bv = BitVector(filename="encrypted_b.txt")

    roundKeys = gen_round_keys(key_bv)
    roundKeys.reverse()

    output = BitVector(size=0)
    # while msg_bv:
    #     if msg_bv.length() < 128:
    #         block = msg_bv.pad_from_right(128 - msg_bv.length())
    #         msg_bv = 0
    #     else:
    #         block, msg_bv = [msg_bv[:128], msg_bv[128:]]
    while msg_bv.more_to_read:
        # fill one block with 128 bits
        block = msg_bv.read_bits_from_file(128)
        # check padding needed
        if block.length() < 128:
            block.pad_from_right(128-block.length())

        # adding 1st round key
        block = Add1stRoundKey(block, roundKeys[0])
        # construct to 4x4 state array
        state_array = [[0 for _ in range(4)] for _ in range(4)]
        for i in range(4):
            for j in range(4):
                state_array[j][i] = BitVector(size=8, intVal=block[32 * i + 8 * j:32 * i + 8 * (j + 1)].int_val())

        # 14 rounds of processing
        for i in range(1, 15):
            # last round
            state_array = AddRoundKey(InvSubBytes(InvShiftRows(state_array)), roundKeys[i])
            if i != 14:
                state_array = InvMixColumns(state_array)


        # convert 4x4 state array to one dimensional array
        for i in range(4):
            for j in range(4):
                output += state_array[i][j]

    with open("decrypted.txt", "w") as file:
        file.write(output.get_text_from_bitvector())


if __name__ == "__main__":
    # encrypt()
    # print("encrypted done")
    decrypt()
