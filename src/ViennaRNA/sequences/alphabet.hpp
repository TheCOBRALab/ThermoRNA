#pragma once

#ifdef VRNA_WARN_DEPRECATED
#if defined(__clang__)
#define DEPRECATED(func, msg) func __attribute__((deprecated("", msg)))
#elif defined(__GNUC__)
#define DEPRECATED(func, msg) func __attribute__((deprecated(msg)))
#else
#define DEPRECATED(func, msg) func
#endif
#else
#define DEPRECATED(func, msg) func
#endif

/**
 *  @file     alphabet.h
 *  @ingroup  utils, alphabet_utils
 *  @brief    Functions to process, convert, and generally handle different nucleotide
 *            and/or base pair alphabets
 */

/**
 *  @addtogroup alphabet_utils
 *  @{
 *  @brief  Functions to cope with various aspects related to the nucleotide sequence alphabet
 */

#include "ViennaRNA/model.hpp"

namespace thermorna::viennarna {

int vrna_nucleotide_IUPAC_identity(char a, char b);

/**
 *  @brief Get a numerical representation of the nucleotide sequence
 *
 *  @param  sequence    The input sequence in upper-case letters
 *  @param  md          A pointer to a #vrna_md_t data structure that specifies the conversion type
 *  @return             A list of integer encodings for each sequence letter (1-based). Position 0
 * denotes the length of the list
 */
short* vrna_seq_encode(const char* sequence, const vrna_md_t* md);

/**
 *  @brief Get a numerical representation of the nucleotide sequence (simple version)
 *
 *  @param  sequence    The input sequence in upper-case letters
 *  @param  md          A pointer to a #vrna_md_t data structure that specifies the conversion type
 *  @return             A list of integer encodings for each sequence letter (1-based). Position 0
 * denotes the length of the list
 */
short* vrna_seq_encode_simple(const char* sequence, const vrna_md_t* md);

/**
 *  @brief  Encode a nucleotide character to numerical value
 *
 *  This function encodes a nucleotide character to its numerical representation as required by many
 * functions in RNAlib.
 *
 *  @see  vrna_nucleotide_decode(), vrna_seq_encode()
 *
 *  @param  c   The nucleotide character to encode
 *  @param  md  The model details that determine the kind of encoding
 *  @return     The encoded nucleotide
 */
int vrna_nucleotide_encode(char c, const vrna_md_t* md);

/**
 *  @brief  Decode a numerical representation of a nucleotide back into nucleotide alphabet
 *
 *  This function decodes a numerical representation of a nucleotide character back into nucleotide
 * alphabet
 *
 *  @see  vrna_nucleotide_encode(), vrna_seq_encode()
 *
 *  @param  enc The encoded nucleotide
 *  @param  md  The model details that determine the kind of decoding
 *  @return     The decoded nucleotide character
 */
char vrna_nucleotide_decode(int enc, const vrna_md_t* md);

void vrna_aln_encode(const char* sequence, short** S_p, short** s5_p, short** s3_p, char** ss_p,
                     unsigned int** as_p, const vrna_md_t* md);

unsigned int vrna_get_ptype_md(int i, int j, const vrna_md_t* md);

unsigned int vrna_get_ptype(int ij, const char* ptype);

unsigned int vrna_get_ptype_window(int i, int j, const char** ptype);

}  // namespace thermorna::viennarna