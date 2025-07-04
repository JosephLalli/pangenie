set (CMAKE_CXX_STANDARD 11)
set (PROGRAM_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)
include_directories (${PROGRAM_SOURCE_DIR})

# Define base project files
file (GLOB_RECURSE  ProjectFiles  ${PROGRAM_SOURCE_DIR}/emissionprobabilitycomputer.cpp ${PROGRAM_SOURCE_DIR}/copynumber.cpp ${PROGRAM_SOURCE_DIR}/kmerpath16.cpp ${PROGRAM_SOURCE_DIR}/kmerpath.cpp ${PROGRAM_SOURCE_DIR}/uniquekmers.cpp ${PROGRAM_SOURCE_DIR}/biallelicuniquekmers.cpp ${PROGRAM_SOURCE_DIR}/multiallelicuniquekmers.cpp ${PROGRAM_SOURCE_DIR}/variant.cpp ${PROGRAM_SOURCE_DIR}/variantreader.cpp ${PROGRAM_SOURCE_DIR}/graphbuilder.cpp ${PROGRAM_SOURCE_DIR}/probabilitycomputer.cpp ${PROGRAM_SOURCE_DIR}/transitionprobabilitycomputer.cpp ${PROGRAM_SOURCE_DIR}/hmm.cpp ${PROGRAM_SOURCE_DIR}/columnindexer.cpp ${PROGRAM_SOURCE_DIR}/columnindexer.cpp ${PROGRAM_SOURCE_DIR}/genotypingresult.cpp ${PROGRAM_SOURCE_DIR}/dnasequence.cpp ${PROGRAM_SOURCE_DIR}/fastareader.cpp ${PROGRAM_SOURCE_DIR}/jellyfishcounter.cpp ${PROGRAM_SOURCE_DIR}/jellyfishreader.cpp ${PROGRAM_SOURCE_DIR}/histogram.cpp ${PROGRAM_SOURCE_DIR}/sequenceutils.cpp ${PROGRAM_SOURCE_DIR}/pathsampler.cpp ${PROGRAM_SOURCE_DIR}/probabilitytable.cpp ${PROGRAM_SOURCE_DIR}/kmerparser.cpp ${PROGRAM_SOURCE_DIR}/graph.cpp ${PROGRAM_SOURCE_DIR}/haplotypesampler.cpp ${PROGRAM_SOURCE_DIR}/samplingemissions.cpp ${PROGRAM_SOURCE_DIR}/samplingtransitions.cpp ${PROGRAM_SOURCE_DIR}/sampledpanel.cpp ${PROGRAM_SOURCE_DIR}/commands.cpp ${PROGRAM_SOURCE_DIR}/commandlineparser.cpp ${PROGRAM_SOURCE_DIR}/timer.cpp ${PROGRAM_SOURCE_DIR}/threadpool.cpp ${PROGRAM_SOURCE_DIR}/stepwiseuniquekmercomputer.cpp ${PROGRAM_SOURCE_DIR}/uniquekmercomputer.cpp ${PROGRAM_SOURCE_DIR}/kmercounter.cpp)

# Add KFF support files if available
if(KFF_FOUND)
    list(APPEND ProjectFiles ${PROGRAM_SOURCE_DIR}/kffreader.cpp)
endif()

# Define test source files
set(TEST_SOURCES 
    tests.cpp utils.cpp EmissionProbabilityComputerTest.cpp CopyNumberTest.cpp 
    UniqueKmersTest.cpp UniqueKmerComputerTest.cpp KmerPathTest.cpp VariantTest.cpp 
    VariantReaderTest.cpp GraphBuilderTest.cpp ProbabilityComputerTest.cpp 
    TransitionProbabilityComputerTest.cpp HMMTest.cpp ColumnIndexerTest.cpp 
    GenotypingResultTest.cpp DnaSequenceTest.cpp FastaReaderTest.cpp KmerCounterTest.cpp 
    HistogramTest.cpp PathSamplerTest.cpp ProbabilityTableTest.cpp KmerParser.cpp 
    HaplotypeSamplerTest.cpp SamplingEmissionsTest.cpp SamplingTransitionsTest.cpp 
    SampledPanelTest.cpp CommandsTest.cpp KffReaderTest.cpp CommandsKffTest.cpp)

add_executable(tests ${TEST_SOURCES} ${ProjectFiles})

# Link libraries
set(TEST_LIBS ${JELLYFISH_LIBRARIES} ${ZLIB_LIBRARIES} ${CEREAL_LIBRARIES})
set(TEST_LDFLAGS ${JELLYFISH_LDFLAGS_OTHER} ${ZLIB_LDFLAGS_OTHER} ${CEREAL_LDFLAGS_OTHER})

if(KFF_FOUND)
    list(APPEND TEST_LIBS ${KFF_LIBRARY})
endif()

target_link_libraries(tests ${TEST_LIBS} ${TEST_LDFLAGS})

# Add KFF test file creation utility if KFF is available
if(KFF_FOUND)
    add_executable(create_minimal_kff create_minimal_kff.cpp)
    target_link_libraries(create_minimal_kff ${KFF_LIBRARY})
    
    # Custom target to create test KFF files
    add_custom_target(create_kff_test_data
        COMMAND ${CMAKE_CURRENT_BINARY_DIR}/create_minimal_kff ${CMAKE_CURRENT_SOURCE_DIR}/data/test.kff
        DEPENDS create_minimal_kff
        COMMENT "Creating test KFF files"
    )
endif()
