execute_process(COMMAND git rev-parse HEAD WORKING_DIRECTORY ${SOURCE} OUTPUT_VARIABLE commit OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND git hash-object src/atom/atom_thunder.h WORKING_DIRECTORY ${SOURCE} OUTPUT_VARIABLE header_blob OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND git hash-object src/atom/atom_thunder.cpp WORKING_DIRECTORY ${SOURCE} OUTPUT_VARIABLE source_blob OUTPUT_STRIP_TRAILING_WHITESPACE)
file(WRITE ${OUTPUT}.next
"#pragma once\n\nnamespace qstorm::version {\ninline constexpr char commit[] = \"${commit}\";\ninline constexpr char header_blob[] = \"${header_blob}\";\ninline constexpr char source_blob[] = \"${source_blob}\";\n}\n")
file(COPY_FILE ${OUTPUT}.next ${OUTPUT} ONLY_IF_DIFFERENT)
