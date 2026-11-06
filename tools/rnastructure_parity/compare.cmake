# Runs both dumpers for every case and checks that ThermoRNA reproduces RNAstructure's tables.
#
#   cmake -DRS_DUMP=<rnastructure_dump> -DTR_DUMP=<thermorna_dump> -DDATA=<data_tables>
#         -DOUT=<scratch dir> [-DCASES="rna:310.15;rna:300;..."] -P compare.cmake
#
# A case is <alphabet>:<temperature in K>. When RNAstructure cannot load a case (e.g. missing
# .dh files), ThermoRNA must refuse it as well.

if(NOT CASES)
    set(CASES
        rna:310.15 dna:310.15 m6A:310.15 PSU:310.15 1:310.15 Full1:310.15 FullP:310.15
        dmem:310.15 msa:310.15
        rna:273.15 rna:298.15 rna:300 rna:323.15 rna:330.5 rna:350 rna:373.15 rna:400
        dna:300 m6A:300 PSU:330)
endif()

file(MAKE_DIRECTORY "${OUT}")
set(failures 0)
foreach(case IN LISTS CASES)
    string(REPLACE ":" ";" parts "${case}")
    list(GET parts 0 alphabet)
    list(GET parts 1 kelvin)
    set(rs_file "${OUT}/rnastructure_${alphabet}_${kelvin}.txt")
    set(tr_file "${OUT}/thermorna_${alphabet}_${kelvin}.txt")
    execute_process(COMMAND "${RS_DUMP}" "${DATA}" "${alphabet}" "${kelvin}"
                    OUTPUT_FILE "${rs_file}" ERROR_QUIET)
    execute_process(COMMAND "${TR_DUMP}" "${DATA}" "${alphabet}" "${kelvin}"
                    OUTPUT_FILE "${tr_file}" ERROR_QUIET)
    file(READ "${rs_file}" rs_text)
    file(READ "${tr_file}" tr_text)
    string(SUBSTRING "${rs_text}" 0 5 rs_head)
    string(SUBSTRING "${tr_text}" 0 5 tr_head)
    if(rs_head STREQUAL "ERROR")
        if(tr_head STREQUAL "ERROR")
            message(STATUS "${case}: refused by both")
        else()
            message(SEND_ERROR "${case}: RNAstructure refuses it but ThermoRNA loads it")
            math(EXPR failures "${failures} + 1")
        endif()
    elseif(rs_text STREQUAL tr_text)
        string(REGEX MATCHALL "\n" lines "${rs_text}")
        list(LENGTH lines count)
        message(STATUS "${case}: identical (${count} lines)")
    else()
        message(SEND_ERROR "${case}: tables differ, compare ${rs_file} and ${tr_file}")
        math(EXPR failures "${failures} + 1")
    endif()
endforeach()

if(failures GREATER 0)
    message(FATAL_ERROR "${failures} case(s) differ")
endif()
