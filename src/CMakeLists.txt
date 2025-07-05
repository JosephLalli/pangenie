set(PANGENIE_SOURCES
	emissionprobabilitycomputer.cpp
	copynumber.cpp
	commandlineparser.cpp
	commands.cpp
	columnindexer.cpp
	dnasequence.cpp
	fastareader.cpp
	genotypingresult.cpp
	graphbuilder.cpp
	graph.cpp
	haplotypesampler.cpp
	histogram.cpp
	hmm.cpp
	jellyfishcounter.cpp
	jellyfishreader.cpp
	kmerpath.cpp
	kmerpath16.cpp
	kmerparser.cpp
	pathsampler.cpp
	probabilitycomputer.cpp
	probabilitytable.cpp
	sampledpanel.cpp
	samplingemissions.cpp
	samplingtransitions.cpp
	sequenceutils.cpp
	stepwiseuniquekmercomputer.cpp
	timer.cpp
	transitionprobabilitycomputer.cpp
	threadpool.cpp
	uniquekmercomputer.cpp
	multiallelicuniquekmers.cpp
	biallelicuniquekmers.cpp
	variant.cpp
	variantreader.cpp)

# Add KFF support if available
if(KFF_FOUND)
	list(APPEND PANGENIE_SOURCES kffreader.cpp)
endif()

add_library(PanGenieLib SHARED ${PANGENIE_SOURCES})


# Define common link libraries
set(COMMON_LIBS PanGenieLib ${JELLYFISH_LIBRARIES} ${ZLIB_LIBRARIES} ${CEREAL_LIBRARIES})
set(COMMON_LDFLAGS ${JELLYFISH_LDFLAGS_OTHER} ${ZLIB_LDFLAGS_OTHER} ${CEREAL_LDFLAGS_OTHER})

if(KFF_FOUND)
    list(APPEND COMMON_LIBS ${KFF_LIBRARY})
endif()

add_executable(PanGenie-index pangenie-index.cpp)
target_link_libraries(PanGenie-index ${COMMON_LIBS} ${COMMON_LDFLAGS})

add_executable(PanGenie pangenie-genotype.cpp)
target_link_libraries(PanGenie ${COMMON_LIBS} ${COMMON_LDFLAGS})

add_executable(PanGenie-vcf pangenie-vcf.cpp)
target_link_libraries(PanGenie-vcf ${COMMON_LIBS} ${COMMON_LDFLAGS})

add_executable(PanGenie-sampling pangenie-sampling.cpp)
target_link_libraries(PanGenie-sampling ${COMMON_LIBS} ${COMMON_LDFLAGS})

add_executable(Analyze-UK analyze-uk.cpp)
target_link_libraries(Analyze-UK ${COMMON_LIBS} ${COMMON_LDFLAGS})
