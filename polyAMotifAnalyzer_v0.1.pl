#!/usr/bin/perl -w

#====================================================================================================================================================#
#<use>
$|++; #---turn on the auto flush for the progress bar
use strict;
use File::Path;
use Time::HiRes qw( time );
use Math::Random;
use Math::Round;
use Storable;
use Getopt::Long;
use File::Basename;
use File::Spec::Functions qw(rel2abs);
use List::Util qw (sum shuffle min max);
use threads;
use threads::shared;
use Statistics::Descriptive;
use URI::Escape;
use Cwd 'abs_path';
#<\use>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<doc>
#	Description
#		This is a Perl script to analyze a predefined set of motifs around the polyA sites.
#
#	Input
#
#		--mRNABasedPolyAInfoHshPlsPath=			path [compulsory]; a hash sotroable contains mRNA Based PolyA Info , mRNABasedPolyAInfoHsh.pls, generated in polyAFinder;
#		--fastaPath=							file path [compulsory]; the path fasta file contains the genome sequence, for generating blank perl storables;
#		--gffPath=								path[compulsory]; path of the reference GFF for gene annotation;
#		--outDir=								output directory
#	
#	v0.1
#		#---[14/10/2013 19:48] 
#
#<\doc>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<lastCmdCalled>
#
#	[2013-10-15 19:23]	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/polyAMotifAnalyzer/v0.1/polyAMotifAnalyzer_v0.1.pl --mRNABasedPolyAInfoHshPlsPath=/Volumes/A_MPro2TB/NGS/analyses/EHI_polyA/polyATailFinderTop3RdALen/IA30.IR20.BR95.BC95.SS200.SA200.MG20.OB500.PB1000.AP1000.MP1000.NP10/storable/mRNABasedPolyAInfoHsh.pls --fastaPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/fasta/genome.sorted.fa --gffPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/gff/forPileupCounter.gff --outDir=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/polyAMotifAnalyzer/v0.1/
#
#	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/polyAMotifAnalyzer/v0.1/polyAMotifAnalyzer_v0.1.pl
#	--mRNABasedPolyAInfoHshPlsPath=/Volumes/A_MPro2TB/NGS/analyses/EHI_polyA/polyATailFinderTop3RdALen/IA30.IR20.BR95.BC95.SS200.SA200.MG20.OB500.PB1000.AP1000.MP1000.NP10/storable/mRNABasedPolyAInfoHsh.pls
#	--fastaPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/fasta/genome.sorted.fa
#	--gffPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/gff/forPileupCounter.gff
#	--outDir=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/polyAMotifAnalyzer/v0.1/
#
#<\lastCmdCalled>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<global>
my $scriptDirPath = dirname(rel2abs($0));
my $globalTmpLogPath = "$scriptDirPath/tmp.log.txt";
open TMPLOG, ">", $globalTmpLogPath;
#<\global>
#====================================================================================================================================================#

#====================================================================================================================================================#
{	#Main sections lexical scope starts
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 0_startingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|2024, readParameters|2288
#	secondaryDependOnSub: currentTime|626
#
#<section ID="startingTasks" num="0">
#-----print start log
&printCMDLogOrFinishMessage("CMDLog");#->2024

#-----read parameters
my ($mRNABasedPolyAInfoHshPlsPath, $fastaPath, $gffPath, $outDir) = &readParameters();#->2288
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 1_defineHardCodedParam
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineHardCodedParam" num="1">
my $siteSeqRng = 150;
my $maxPolymerSize = 5;
my $abvBkgrdFactor = 1.2; #----the background will be scaled up using this factor during comparison. The purpose is to indentify the signal thay is truely above background
my $extraMargin = 0; #---extra margin added to the motif boundarys defined
my $logTransformPVal = "yes"; #-----use 1/pValue or log(1/pValue);
my $minStreakHsh_ref = {}; #---minimum nt streak to start and end the motif boundary 
$minStreakHsh_ref->{'positive'} = 3;
$minStreakHsh_ref->{'negative'} = 3;
my $minCDSLength = 1000; #---minimum CDS length to be plot in cov along CDS
my $maxThread = 8;
my $PAScoreAlongCDSBinWidth = 5;

my $predictMotifInfoHsh_ref = {};
$predictMotifInfoHsh_ref->{'CSE'}{'scoreFactor'} = 1; #---the factor use to scale the score of the motif
$predictMotifInfoHsh_ref->{'PAS'}{'scoreFactor'} = 1; #---the factor use to scale the score of the motif
$predictMotifInfoHsh_ref->{'DSE'}{'scoreFactor'} = 1; #---the factor use to scale the score of the motif

$predictMotifInfoHsh_ref->{'CSE'}{'maxPValGenomePredict'} = 0.005; #---the max pval of motif to be valid
$predictMotifInfoHsh_ref->{'PAS'}{'maxPValGenomePredict'} = 0.05; #---the max pval of motif to be valid
$predictMotifInfoHsh_ref->{'DSE'}{'maxPValGenomePredict'} = 0.001; #---the max pval of motif to be valid

$predictMotifInfoHsh_ref->{'CSE'}{'maxPValDefineBound'} = 0.005; #---the max pval of motif to be valid for defining the positional bound
$predictMotifInfoHsh_ref->{'PAS'}{'maxPValDefineBound'} = 0.005; #---the max pval of motif to be valid for defining the positional bound
$predictMotifInfoHsh_ref->{'DSE'}{'maxPValDefineBound'} = 0.005; #---the max pval of motif to be valid for defining the positional bound

$predictMotifInfoHsh_ref->{'CSE'}{'mustValid'} = 'yes'; #---the motif must be valid or not
$predictMotifInfoHsh_ref->{'PAS'}{'mustValid'} = 'yes'; #---the motif must be valid or not
$predictMotifInfoHsh_ref->{'DSE'}{'mustValid'} = 'yes'; #---the motif must be valid or not
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 2_defineHardCodedDir
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineHardCodedDir" num="2">
my @resultDirTagAry = ();
foreach my $motif (qw/PAS CSE DSE/) {
	push @resultDirTagAry, "$motif.$predictMotifInfoHsh_ref->{$motif}{'maxPValDefineBound'}";
}

my %PAScoreResultTagHsh = ();
foreach my $motif (keys %{$predictMotifInfoHsh_ref}) {
	my @tmpAry = ();
	foreach my $param (qw/scoreFactor maxPValGenomePredict mustValid/) {
		push @tmpAry, $predictMotifInfoHsh_ref->{$motif}{$param};
	}
	$PAScoreResultTagHsh{$motif} = join "_", @tmpAry;
}
	
my $PAScoreResultTag = "PAS.$PAScoreResultTagHsh{'PAS'}.CSE.$PAScoreResultTagHsh{'CSE'}.DSE.$PAScoreResultTagHsh{'DSE'}.logPVal.$logTransformPVal";

my $resultDirTag = join ".", @resultDirTagAry;
my @mkDirAry;
my $PAScoreDir = "$outDir/$resultDirTag/PAScore/$PAScoreResultTag/"; push @mkDirAry, $PAScoreDir;
my $PAScoreStorableDir = "$PAScoreDir/storable"; push @mkDirAry, $PAScoreStorableDir;
my $PAScoreWigDir = "$PAScoreDir/wig"; push @mkDirAry, $PAScoreWigDir;
my $PAScoreLogDir = "$PAScoreDir/log"; push @mkDirAry, $PAScoreLogDir;
my $resultStorableDir = "$outDir/$resultDirTag/storable/"; push @mkDirAry, $resultStorableDir;
my $resultLogDir = "$outDir/$resultDirTag/log/"; push @mkDirAry, $resultLogDir;
my $resultMotifDir = "$outDir/$resultDirTag/motifUsed/"; push @mkDirAry, $resultMotifDir;
my $resultGFFDir = "$outDir/$resultDirTag/GFF/"; push @mkDirAry, $resultGFFDir;
my $resultWigDir = "$outDir/$resultDirTag/wig/"; push @mkDirAry, $resultWigDir;
my $resultXMLDir = "$outDir/$resultDirTag/XML/"; push @mkDirAry, $resultXMLDir;
my $resultFastaDir = "$outDir/$resultDirTag/fasta/"; push @mkDirAry, $resultFastaDir;
my $mastDir = "$outDir/$resultDirTag/MAST/"; push @mkDirAry, $mastDir;
my $mastBackgroundDir = "$mastDir/background/"; push @mkDirAry, $mastBackgroundDir;
my $mastMotifDir = "$mastDir/motif/"; push @mkDirAry, $mastMotifDir;
my $mastRunDir = "$mastDir/run/"; push @mkDirAry, $mastRunDir;

my $generalggplotDirHsh_ref = {};
my $PAScoreggplotDirHsh_ref = {};
foreach my $fileType (qw /dat pdf R log/) {
	$generalggplotDirHsh_ref->{$fileType} = "$outDir/$resultDirTag/ggplot/$fileType"; push @mkDirAry, $generalggplotDirHsh_ref->{$fileType};
	$PAScoreggplotDirHsh_ref->{$fileType} = "$PAScoreDir/ggplot/$fileType"; push @mkDirAry, $PAScoreggplotDirHsh_ref->{$fileType};
}

my $weblogoDirHsh_ref = {};
foreach my $fileType (qw /pdf fasta cmd/) {$weblogoDirHsh_ref->{$fileType} = "$outDir/$resultDirTag/weblogo/$fileType"; push @mkDirAry, $weblogoDirHsh_ref->{$fileType};}

system ("mkdir -pm 777 $_") foreach @mkDirAry;
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 3_defineOutFilePath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineOutFilePath" num="3">
#---[14/10/2013 17:26] the motifs are from polyAFinder
my $motifFilePathHsh_ref = {};
$motifFilePathHsh_ref->{'CSE'} = "/Volumes/A_MPro2TB/NGS/analyses/EHI_polyA/polyATailFinderTop3RdALen/IA30.IR20.BR95.BC95.SS200.SA200.MG20.OB500.PB1000.AP1000.MP1000.NP10/motif/s_tail.CSE.motif.txt";
$motifFilePathHsh_ref->{'DSE'} = "/Volumes/A_MPro2TB/NGS/analyses/EHI_polyA/polyATailFinderTop3RdALen/IA30.IR20.BR95.BC95.SS200.SA200.MG20.OB500.PB1000.AP1000.MP1000.NP10/motif/s_tail.DSE.negative.motif.txt";
$motifFilePathHsh_ref->{'PAS'} = "/Volumes/A_MPro2TB/NGS/analyses/EHI_polyA/polyATailFinderTop3RdALen/IA30.IR20.BR95.BC95.SS200.SA200.MG20.OB500.PB1000.AP1000.MP1000.NP10/motif/s_tail.PAS.negative.motif.txt";

my $bkgdNtFreqHshPath = "$resultStorableDir/bkgdNtFreqHsh.pls";

#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 4_processInputData
#	primaryDependOnSub: checkGeneInfo|464, generateFastaLengthHsh|699, getCtgryGeneInfo|957, readGFF_oneRNAPerGene|2131, readMultiFasta|2230, reverseComplementRefFasta|2343
#	secondaryDependOnSub: currentTime|626, reportStatus|2322
#
#<section ID="processInputData" num="4">

#----------Read fasta
my $fastaHsh_ref = &readMultiFasta($fastaPath);#->2230
my ($fastaWithRevComHsh_ref, $revComFastaPath) = &reverseComplementRefFasta($fastaHsh_ref, $resultFastaDir);#->2343
my ($fastaLengthHsh_ref) = &generateFastaLengthHsh($fastaHsh_ref);#->699

#----------Read Gff
my ($geneInfoHsh_ref) = &readGFF_oneRNAPerGene($gffPath);#->2131
&checkGeneInfo($geneInfoHsh_ref);#->464

#----------Get bfmRNA and mRNA ranges
my @mRNAAry = qw/bfmRNA/;
my ($mRNAInfoHsh_ref, $mRNAByCntgHsh_ref)= &getCtgryGeneInfo($geneInfoHsh_ref, \@mRNAAry);#->957
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 5_retrieveSequenceSurroundingPredefinedTSS
#	primaryDependOnSub: calculateBackgroundNucleotideFrequency|374, generateShuffleSeq|763, getSequenceAroundmRNAReferencePoint|1267, getmRNAReferencePoints|1374, plotBaseCompositionAroundmRNAReferencePoint|1572, plotWeblogoAroundRefPoint|1863
#	secondaryDependOnSub: calculateBaseCompositionInAlignments|420, createWeblogo|600, generateMASTBackgroundFile|718, ggplotXYLinesMultipleSamples|1537, reportStatus|2322
#
#<section ID="retrieveSequenceSurroundingPredefinedTSS" num="5">
my $mRNABasedPolyAInfoHsh_ref = retrieve($mRNABasedPolyAInfoHshPlsPath);

my ($mRNARefPtHsh_ref) = &getmRNAReferencePoints($mRNAInfoHsh_ref, $mRNABasedPolyAInfoHsh_ref);#->1374

my ($seqAroundSiteHsh_ref, $seqAroundSiteInfoHsh_ref) = &getSequenceAroundmRNAReferencePoint($mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $fastaHsh_ref, $siteSeqRng, $resultFastaDir, $resultStorableDir);#->1267

&plotBaseCompositionAroundmRNAReferencePoint($seqAroundSiteHsh_ref, $generalggplotDirHsh_ref);#->1572
&plotWeblogoAroundRefPoint($seqAroundSiteHsh_ref, $weblogoDirHsh_ref);#->1863

my ($shuffleSeqAroundSiteHsh_ref, $shuffleSeqAroundSiteInfoHsh_ref) = &generateShuffleSeq($seqAroundSiteHsh_ref, $resultFastaDir, $resultStorableDir);#->763

my ($bkgdNtFreqHsh_ref) = &calculateBackgroundNucleotideFrequency($seqAroundSiteHsh_ref, $maxPolymerSize, $mastBackgroundDir, $bkgdNtFreqHshPath);#->374

#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 6_scanMotifOccurence
#	primaryDependOnSub: scanMotifAroundSiteWithMAST|2395, scanMotifWholeGenomeWithMAST|2454
#	secondaryDependOnSub: createMotifHitHshStorable|561, generateMASTBackgroundFile|718, getMastGenomeBothStrandHit|994, getSingleMotifMASTLogPostionalData|1337, reportStatus|2322, runMAST|2374, storeMotifHitToHshStorable|2501
#
#<section ID="scanMotifOccurence" num="6">
my ($cntgMotifHitPlsPathHsh_ref) = &scanMotifWholeGenomeWithMAST($revComFastaPath, $fastaHsh_ref, $motifFilePathHsh_ref, $mastRunDir, $maxPolymerSize, $mastBackgroundDir, $fastaLengthHsh_ref, $resultStorableDir);#->2454

my ($mastAroundSiteResultHsh_ref) = &scanMotifAroundSiteWithMAST($seqAroundSiteInfoHsh_ref, $shuffleSeqAroundSiteInfoHsh_ref, $motifFilePathHsh_ref, $bkgdNtFreqHsh_ref, $mastRunDir, $predictMotifInfoHsh_ref, $resultStorableDir);#->2395

#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 7_defineMotifBoundary
#	primaryDependOnSub: defineMotifPositionBoundsUsingShuffleBackground|644, plotMotifPostionFactor|1724
#	secondaryDependOnSub: getMotifBoundCutoff|1031, ggplotXYLinesMultipleSamples|1537
#
#<section ID="defineMotifBoundary" num="7">
my ($motifPostionFactorHsh_ref) = &defineMotifPositionBoundsUsingShuffleBackground($predictMotifInfoHsh_ref, $mastAroundSiteResultHsh_ref, $generalggplotDirHsh_ref, $minStreakHsh_ref, $abvBkgrdFactor, $extraMargin, $siteSeqRng);#->644
&plotMotifPostionFactor($motifPostionFactorHsh_ref, $generalggplotDirHsh_ref);#->1724

#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 8_destroyBigUnusedVar
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="destroyBigUnusedVar" num="8">
undef $mastAroundSiteResultHsh_ref;
undef $seqAroundSiteHsh_ref;
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 9_predictTSS
#	primaryDependOnSub: predictGenomeWidePolyASite|1903, printPAScoreWiggle|2057
#	secondaryDependOnSub: checkRunningThreadAndWaitToJoin|490, createEmptyStorableForGenowideTSSPredictionData|518, generateThreadHshWithRandomCntg|825, printWiggleSingleTrackFromCntgCovPlsPathHsh|2093, reportStatus|2322
#
#<section ID="predictTSS" num="9">
my ($genomeWidePolyASitePlsPathHsh_ref) = &predictGenomeWidePolyASite($motifPostionFactorHsh_ref, $cntgMotifHitPlsPathHsh_ref, $predictMotifInfoHsh_ref, $fastaHsh_ref, $PAScoreStorableDir, $logTransformPVal);#->1903
&printPAScoreWiggle($PAScoreWigDir, $genomeWidePolyASitePlsPathHsh_ref);#->2057
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 10_analyzePAScore
#	primaryDependOnSub: getPAScoreAroundmRNAReferencePoint|1113, getPAScoreInExonAndTSS|1207, plotCovAlongCDSBothStrnd|1628, plotPAScoreAroundmRNAReferencePoint|1761, plotPAScoreInExonAndTSS|1805
#	secondaryDependOnSub: getCoverageOfItemRngType_multiStrand|852, ggplotMultiSampleBoxWhisker|1418, ggplotTwoSampleHistogram|1463, ggplotXYLinesMultipleSamples|1537, reportStatus|2322
#
#<section ID="analyzePAScore" num="10">
my ($PAScoreExonTSSHshNonZero_ref, $PAScoreExonTSSHshWithZero_ref) = &getPAScoreInExonAndTSS($genomeWidePolyASitePlsPathHsh_ref, $mRNAByCntgHsh_ref, $mRNARefPtHsh_ref, $siteSeqRng, $mRNAInfoHsh_ref);#->1207
&plotPAScoreInExonAndTSS($PAScoreExonTSSHshNonZero_ref, $PAScoreExonTSSHshWithZero_ref, $PAScoreggplotDirHsh_ref);#->1805

my ($PAScoreRefPtPlotHsh_ref, $PAScoreRefPtBymRNAHsh_ref) = &getPAScoreAroundmRNAReferencePoint($mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $genomeWidePolyASitePlsPathHsh_ref, $siteSeqRng, $mRNAByCntgHsh_ref);#->1113

&plotPAScoreAroundmRNAReferencePoint($PAScoreRefPtPlotHsh_ref, $PAScoreggplotDirHsh_ref, $mRNARefPtHsh_ref);#->1761

(undef, undef) = &plotCovAlongCDSBothStrnd($genomeWidePolyASitePlsPathHsh_ref, $mRNAInfoHsh_ref, $mRNAByCntgHsh_ref, $minCDSLength, $PAScoreggplotDirHsh_ref, $maxThread, $PAScoreAlongCDSBinWidth);#->1628

#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 11_finishingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|2024
#	secondaryDependOnSub: currentTime|626
#
#<section ID="finishingTasks" num="11">
&printCMDLogOrFinishMessage("finishMessage");#->2024
close TMPLOG;
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
}	#Main sections lexical scope ends
#====================================================================================================================================================#

#====================================================================================================================================================#
#List of subroutines by category
#
#	alignment [n=1]:
#		calculateBaseCompositionInAlignments
#
#	baseComposition [n=1]:
#		calculateBaseCompositionInAlignments
#
#	coverage [n=2]:
#		getCoverageOfItemRngType_multiStrand, plotCovAlongCDSBothStrnd
#
#	fasta [n=1]:
#		readMultiFasta
#
#	general [n=8]:
#		checkGeneInfo, currentTime, getCtgryGeneInfo
#		printCMDLogOrFinishMessage, readGFF_oneRNAPerGene, readMultiFasta
#		readParameters, reportStatus
#
#	gff [n=3]:
#		checkGeneInfo, getCtgryGeneInfo, readGFF_oneRNAPerGene
#
#	ggplot [n=2]:
#		ggplotTwoSampleHistogram, ggplotXYLinesMultipleSamples
#
#	multithread [n=2]:
#		checkRunningThreadAndWaitToJoin, generateThreadHshWithRandomCntg
#
#	plotInR [n=2]:
#		ggplotTwoSampleHistogram, ggplotXYLinesMultipleSamples
#
#	reporting [n=1]:
#		currentTime
#
#	specific [n=2]:
#		getmRNAReferencePoints, plotWeblogoAroundRefPoint
#
#	thirdPartyApp [n=1]:
#		plotWeblogoAroundRefPoint
#
#	thridPartyApp [n=2]:
#		createWeblogo, runMAST
#
#	unassigned [n=24]:
#		calculateBackgroundNucleotideFrequency, createEmptyStorableForGenowideTSSPredictionData, createMotifHitHshStorable
#		defineMotifPositionBoundsUsingShuffleBackground, generateFastaLengthHsh, generateMASTBackgroundFile
#		generateShuffleSeq, getMastGenomeBothStrandHit, getMotifBoundCutoff
#		getPAScoreAroundmRNAReferencePoint, getPAScoreInExonAndTSS, getSequenceAroundmRNAReferencePoint
#		getSingleMotifMASTLogPostionalData, ggplotMultiSampleBoxWhisker, plotBaseCompositionAroundmRNAReferencePoint
#		plotMotifPostionFactor, plotPAScoreAroundmRNAReferencePoint, plotPAScoreInExonAndTSS
#		predictGenomeWidePolyASite, printPAScoreWiggle, reverseComplementRefFasta
#		scanMotifAroundSiteWithMAST, scanMotifWholeGenomeWithMAST, storeMotifHitToHshStorable
#
#	wiggle [n=1]:
#		printWiggleSingleTrackFromCntgCovPlsPathHsh
#
#====================================================================================================================================================#

sub calculateBackgroundNucleotideFrequency {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: generateMASTBackgroundFile|718, reportStatus|2322
#	appearInSub: >none
#	primaryAppearInSection: 5_retrieveSequenceSurroundingPredefinedTSS|213
#	secondaryAppearInSection: >none
#	input: $bkgdNtFreqHshPath, $mastBackgroundDir, $maxPolymerSize, $seqAroundSiteHsh_ref
#	output: $bkgdNtFreqHsh_ref
#	toCall: my ($bkgdNtFreqHsh_ref) = &calculateBackgroundNucleotideFrequency($seqAroundSiteHsh_ref, $maxPolymerSize, $mastBackgroundDir, $bkgdNtFreqHshPath);
#	calledInLine: 229
#....................................................................................................................................................#
	my ($seqAroundSiteHsh_ref, $maxPolymerSize, $mastBackgroundDir, $bkgdNtFreqHshPath) = @_;
	
	my $bkgdNtFreqHsh_ref = {};
	
	if (-s $bkgdNtFreqHshPath) {
	
		$bkgdNtFreqHsh_ref = retrieve($bkgdNtFreqHshPath);
		&reportStatus("Retrieving nucleotide frequencies around all siteTypes", 0,"\n");#->2322
	
	} else {
		foreach my $siteType (keys %{$seqAroundSiteHsh_ref}) {
			my $siteTypeFullSeqHsh_ref = {};
			&reportStatus("Getting nucleotide frequencies around $siteType", 0,"\n");#->2322
		
			foreach my $mRNAID (keys %{$seqAroundSiteHsh_ref->{$siteType}}) {
				$siteTypeFullSeqHsh_ref->{'full'}{$mRNAID} = $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'}.$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'};
				$siteTypeFullSeqHsh_ref->{'upStrm'}{$mRNAID} = $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'};
				$siteTypeFullSeqHsh_ref->{'dnStrm'}{$mRNAID} = $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'};
			}

			foreach my $fullOrUpStrmOrDnStrm (keys %{$siteTypeFullSeqHsh_ref}) {
				my $seqHsh_ref = \%{$siteTypeFullSeqHsh_ref->{$fullOrUpStrmOrDnStrm}};
				my $bfilePath = "$mastBackgroundDir/$siteType.$fullOrUpStrmOrDnStrm.mast.bfile.txt";
				my ($freqHsh_ref, $proportionHsh_ref) = &generateMASTBackgroundFile($seqHsh_ref, $bfilePath, $maxPolymerSize);#->718
				$bkgdNtFreqHsh_ref->{$siteType}{$fullOrUpStrmOrDnStrm}{'proportionHsh_ref'} = $proportionHsh_ref;
				$bkgdNtFreqHsh_ref->{$siteType}{$fullOrUpStrmOrDnStrm}{'freqHsh_ref'} = $freqHsh_ref;
				$bkgdNtFreqHsh_ref->{$siteType}{$fullOrUpStrmOrDnStrm}{'bfilePath'} = $bfilePath;
			}
		}
		store($bkgdNtFreqHsh_ref, $bkgdNtFreqHshPath);
	}

	return ($bkgdNtFreqHsh_ref);
}
sub calculateBaseCompositionInAlignments {
#....................................................................................................................................................#
#	subroutineCategory: alignment, baseComposition
#	dependOnSub: reportStatus|2322
#	appearInSub: plotBaseCompositionAroundmRNAReferencePoint|1572
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 5_retrieveSequenceSurroundingPredefinedTSS|213
#	input: $seqAlignHsh_ref
#	output: $baseCompByBaseHsh_ref
#	toCall: my ($baseCompByBaseHsh_ref) = &calculateBaseCompositionInAlignments($seqAlignHsh_ref);
#	calledInLine: 1609
#....................................................................................................................................................#

	my ($seqAlignHsh_ref) = @_;
	
	my $baseCountByPosHsh_ref = {};
	my $baseCompByBaseHsh_ref = {};
	my $validSeqNum = 0;
	my $tmpLengthHsh_ref = {};
	
	foreach my $seqName (keys %{$seqAlignHsh_ref}) {
		next if $seqAlignHsh_ref->{$seqName} =~ m/[^ATGCatgc]/;
		$validSeqNum++;
		my @seqAry = split //, $seqAlignHsh_ref->{$seqName};
		$tmpLengthHsh_ref->{@seqAry}++;
		for my $pos (0..$#seqAry) {
			my $base = $seqAry[$pos];
			$base =~ tr/atgc/ATGC/;
			$baseCountByPosHsh_ref->{$pos}{$base}++;
		}
	}
	
	my $lengthNum = keys %{$tmpLengthHsh_ref};
	&reportStatus("WARNING: Length of the sequences in the alignment is not uniform", 10, "\n") if $lengthNum > 1;#->2322
	foreach my $base (qw/A T G C/) {
		foreach my $pos (sort {$a <=> $b} keys %{$baseCountByPosHsh_ref}) {
			$baseCountByPosHsh_ref->{$pos}{$base} = 0 if not $baseCountByPosHsh_ref->{$pos}{$base};
			$baseCompByBaseHsh_ref->{$base}{$pos} = $baseCountByPosHsh_ref->{$pos}{$base}/$validSeqNum;
		}
	}
	
	return $baseCompByBaseHsh_ref;
	
}
sub checkGeneInfo {
#....................................................................................................................................................#
#	subroutineCategory: general, gff
#	dependOnSub: reportStatus|2322
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|191
#	secondaryAppearInSection: >none
#	input: $geneInfoHsh_ref
#	output: none
#	toCall: &checkGeneInfo($geneInfoHsh_ref);
#	calledInLine: 204
#....................................................................................................................................................#
	
	my ($geneInfoHsh_ref) = @_;
	
	&reportStatus("Checking gene categories", 0, "\n");#->2322
	my $ctrgyCountHsh_ref = {};
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {
		my $ctgry = $geneInfoHsh_ref->{$geneID}{'ctgry'};
		$ctrgyCountHsh_ref->{$ctgry}++;
	}
	
	foreach my $ctgry (sort keys %{$ctrgyCountHsh_ref}) {
		&reportStatus("Item in $ctgry = $ctrgyCountHsh_ref->{$ctgry}", 0, "\n");#->2322
	}
}
sub checkRunningThreadAndWaitToJoin {
#....................................................................................................................................................#
#	subroutineCategory: multithread
#	dependOnSub: reportStatus|2322
#	appearInSub: predictGenomeWidePolyASite|1903, printPAScoreWiggle|2057
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 9_predictTSS|271
#	input: $sleepTime, $verbose
#	output: none
#	toCall: &checkRunningThreadAndWaitToJoin($verbose, $sleepTime);
#	calledInLine: 2013, 2088
#....................................................................................................................................................#
	
	my ($verbose, $sleepTime) = @_;
	
	my @runningThrAry = threads->list(threads::running);
	my @joinableThrAry = threads->list(threads::joinable);
	while (@runningThrAry or @joinableThrAry) {
		@runningThrAry = threads->list(threads::running);
		@joinableThrAry = threads->list(threads::joinable);
		foreach my $joinableThr (@joinableThrAry) {
			$joinableThr->detach() if not $joinableThr->is_running();
		}
		my $numThreadRunning = scalar @runningThrAry;
		&reportStatus("The last $numThreadRunning threads are still running", 20, "\r") if $verbose eq 'yes';#->2322
		sleep $sleepTime;
	}
}
sub createEmptyStorableForGenowideTSSPredictionData {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|2322
#	appearInSub: predictGenomeWidePolyASite|1903
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 9_predictTSS|271
#	input: $fastaHsh_ref, $resultStorableDir
#	output: $allGenomeWideTSSPlsPathExist, $genomeWidePolyASitePlsPathHsh_ref
#	toCall: my ($genomeWidePolyASitePlsPathHsh_ref, $allGenomeWideTSSPlsPathExist) = &createEmptyStorableForGenowideTSSPredictionData($fastaHsh_ref, $resultStorableDir);
#	calledInLine: 1916
#....................................................................................................................................................#
	my ($fastaHsh_ref, $resultStorableDir) = @_;
	
	my $genomeWidePolyASitePlsDir = "$resultStorableDir/genomeWidePolyASitePrediction/";
	system ("mkdir -pm 777 $genomeWidePolyASitePlsDir");

	my $genomeWidePolyASitePlsIdxHsh_ref = {};
	my $genomeWidePolyASitePlsPathHsh_ref = {};
	my $allGenomeWideTSSPlsPathExist = 'yes';
	
	foreach my $cntg (keys %{$fastaHsh_ref}) {
		my $cntgPlsName = "$cntg.ary.pls";
		my $cntgPlsPath = "$genomeWidePolyASitePlsDir/$cntgPlsName";
		$genomeWidePolyASitePlsIdxHsh_ref->{$cntg} = $cntgPlsName;
		$genomeWidePolyASitePlsPathHsh_ref->{$cntg} = $cntgPlsPath;
		if (not -s $cntgPlsPath) {

			&reportStatus("Creating empty storabe for $cntg", 20, "\r");#->2322

			my $cntgPosAry_ref = ();
			push @{$cntgPosAry_ref}, undef foreach (1..length($fastaHsh_ref->{$cntg}));

			store($cntgPosAry_ref, $cntgPlsPath);
			$allGenomeWideTSSPlsPathExist = 'no';
		}
	}

	my $genomeWidePolyASitePlsIdxPath = "$genomeWidePolyASitePlsDir/index.hsh.pls";
	store($genomeWidePolyASitePlsIdxHsh_ref, $genomeWidePolyASitePlsIdxPath);

	return ($genomeWidePolyASitePlsPathHsh_ref, $allGenomeWideTSSPlsPathExist);
}
sub createMotifHitHshStorable {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|2322
#	appearInSub: scanMotifWholeGenomeWithMAST|2454
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_scanMotifOccurence|235
#	input: $fastaHsh_ref, $resultStorableDir
#	output: $allCntgMotifHitPlsExist, $cntgMotifHitIdxPlsPath, $cntgMotifHitPlsPathHsh_ref
#	toCall: my ($cntgMotifHitIdxPlsPath, $cntgMotifHitPlsPathHsh_ref, $allCntgMotifHitPlsExist) = &createMotifHitHshStorable($fastaHsh_ref, $resultStorableDir);
#	calledInLine: 2467
#....................................................................................................................................................#
	my ($fastaHsh_ref, $resultStorableDir) = @_;
	
	my $cntgMotifHitHshPlsDir = "$resultStorableDir/motifHitHsh/";
	my $cntgMotifHitIdxPlsPath = "$cntgMotifHitHshPlsDir/index.hsh.pls";
	my $allCntgMotifHitPlsExist = 'yes';
	
	&reportStatus("Creating empty motif hit hsh storable for motifHits", 0,"\n");#->2322

	system ("mkdir -pm 777 $cntgMotifHitHshPlsDir");
	my $cntgMotifHitPlsIdxHsh_ref = {};
	my $cntgMotifHitPlsPathHsh_ref = {};
	foreach my $cntg (keys %{$fastaHsh_ref}) {
		my $cntgMotifHitPlsName = "$cntg.motifHitHsh.pls";
		my $cntgMotifHitPlsPath = "$cntgMotifHitHshPlsDir/$cntgMotifHitPlsName";
		$cntgMotifHitPlsPathHsh_ref->{$cntg} = $cntgMotifHitPlsPath;
		$cntgMotifHitPlsIdxHsh_ref->{$cntg} = $cntgMotifHitPlsName;

		if (not -s $cntgMotifHitPlsPath) {
			my $cntgMotifHitHsh_ref = {};
			store($cntgMotifHitHsh_ref, $cntgMotifHitPlsPath);
			$allCntgMotifHitPlsExist = 'no';
		}
	}
	store($cntgMotifHitPlsIdxHsh_ref, $cntgMotifHitIdxPlsPath);
	
	return ($cntgMotifHitIdxPlsPath, $cntgMotifHitPlsPathHsh_ref, $allCntgMotifHitPlsExist);
}
sub createWeblogo {
#....................................................................................................................................................#
#	subroutineCategory: thridPartyApp
#	dependOnSub: >none
#	appearInSub: plotWeblogoAroundRefPoint|1863
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 5_retrieveSequenceSurroundingPredefinedTSS|213
#	input: $fastaPath, $pdfPath, $seqAlignHsh_ref, $seqType, $title
#	output: none
#	toCall: &createWeblogo($seqAlignHsh_ref, $pdfPath, $fastaPath, $seqType, $title);
#	calledInLine: 1897
#....................................................................................................................................................#

	my ($seqAlignHsh_ref, $pdfPath, $fastaPath, $seqType, $title) = @_;
	
	#---print the fasta
	open (FASTA, ">", $fastaPath);
	foreach my $seqName (keys %{$seqAlignHsh_ref}) {
		print FASTA ">$seqName\n";
		print FASTA "$seqAlignHsh_ref->{$seqName}\n";
	}
	close FASTA;
	
	system ("weblogo --fin $fastaPath --datatype fasta --format pdf --fout $pdfPath --sequence-type $seqType --title \"$title\" --color-scheme classic");
	
}
sub currentTime {
#....................................................................................................................................................#
#	subroutineCategory: general, reporting
#	dependOnSub: >none
#	appearInSub: printCMDLogOrFinishMessage|2024, readGFF_oneRNAPerGene|2131, reportStatus|2322
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 0_startingTasks|69, 11_finishingTasks|300, 4_processInputData|191
#	input: none
#	output: $runTime
#	toCall: my ($runTime) = &currentTime();
#	calledInLine: 2044, 2047, 2052, 2151, 2338
#....................................................................................................................................................#
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
	
	return $runTime;
}
sub defineMotifPositionBoundsUsingShuffleBackground {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: getMotifBoundCutoff|1031, ggplotXYLinesMultipleSamples|1537
#	appearInSub: >none
#	primaryAppearInSection: 7_defineMotifBoundary|248
#	secondaryAppearInSection: >none
#	input: $abvBkgrdFactor, $extraMargin, $generalggplotDirHsh_ref, $mastAroundSiteResultHsh_ref, $minStreakHsh_ref, $predictMotifInfoHsh_ref, $siteSeqRng
#	output: $motifPostionFactorHsh_ref
#	toCall: my ($motifPostionFactorHsh_ref) = &defineMotifPositionBoundsUsingShuffleBackground($predictMotifInfoHsh_ref, $mastAroundSiteResultHsh_ref, $generalggplotDirHsh_ref, $minStreakHsh_ref, $abvBkgrdFactor, $extraMargin, $siteSeqRng);
#	calledInLine: 253
#....................................................................................................................................................#
	my ($predictMotifInfoHsh_ref, $mastAroundSiteResultHsh_ref, $generalggplotDirHsh_ref, $minStreakHsh_ref, $abvBkgrdFactor, $extraMargin, $siteSeqRng) = @_;

	my $motifPostionFactorHsh_ref = {};
	foreach my $siteType ('polyA_tail') {
		foreach my $motif (keys %{$mastAroundSiteResultHsh_ref->{$siteType}}) {
			my $bothMotifPctHsh_ref = {};
			my $posPctHsh_ref = {};
			foreach my $queryOrShuffle (keys %{$mastAroundSiteResultHsh_ref->{$siteType}{$motif}}) {
				my $motifPctHsh_ref = $mastAroundSiteResultHsh_ref->{$siteType}{$motif}{$queryOrShuffle}{'motifPctHsh_ref'};#---take the has out of the lexcial scope
				foreach my $pos (keys %{$motifPctHsh_ref}) {
					$bothMotifPctHsh_ref->{$queryOrShuffle}{$pos} = $motifPctHsh_ref->{$pos};
					$posPctHsh_ref->{$pos}{$queryOrShuffle} = $motifPctHsh_ref->{$pos};
				}
			}

			my ($leftBound, $rightBound, $postionFactorHsh_ref) = &getMotifBoundCutoff($posPctHsh_ref, $minStreakHsh_ref, $abvBkgrdFactor, $extraMargin, $siteSeqRng, $motif);#->1031
			$motifPostionFactorHsh_ref->{$motif} = $postionFactorHsh_ref;
			
			#print $siteType."\t".$motif."\t".$leftBound."\t".$rightBound."\n";
			my @posAry = (sort {$a <=> $b} keys %{$posPctHsh_ref});
			my $maxPos = $posAry[-1];
			my $maxHitPVal = $predictMotifInfoHsh_ref->{$motif}{'maxPValDefineBound'};
			my $nameTag = "mast.$motif.p$maxHitPVal.in.$siteType";
			my $plotDataHsh_ref = $bothMotifPctHsh_ref;
			my $pdfPath = $generalggplotDirHsh_ref->{'pdf'}."/$nameTag.pdf";
			my $dataPath = $generalggplotDirHsh_ref->{'dat'}."/$nameTag.dat";
			my $RScriptPath = $generalggplotDirHsh_ref->{'R'}."/$nameTag.R";
			my $logPath = $generalggplotDirHsh_ref->{'log'}."/$nameTag.log";
			my $xAxis = 'relativePositon';
			my $YAxis = 'proportion';
			my $YVariable = 'motif';
			#my $extraArg = '+ ylim(0, 100)';
			my $extraArg = " + scale_x_continuous(breaks=seq(0, $maxPos, by=10))";
			$extraArg .= " + geom_vline(xintercept=c($leftBound), linetype=\"dotted\") + annotate(\"text\", x=$leftBound, y=0, label=\"$leftBound\", vjust=-0.2, hjust=-0.1, angle=90)";
			$extraArg .= " + geom_vline(xintercept=c($rightBound), linetype=\"dotted\") + annotate(\"text\", x=$rightBound, y=0, label=\"$rightBound\", vjust=-0.2, hjust=-0.1, angle=90)";
			my $height = 6;
			my $width = 14;
			&ggplotXYLinesMultipleSamples($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $xAxis, $YAxis, $YVariable, $extraArg, $height, $width);#->1537
		}
	}

	return ($motifPostionFactorHsh_ref);
}
sub generateFastaLengthHsh {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|191
#	secondaryAppearInSection: >none
#	input: $fastaHsh_ref
#	output: $fastaLengthHsh_ref
#	toCall: my ($fastaLengthHsh_ref) = &generateFastaLengthHsh($fastaHsh_ref);
#	calledInLine: 200
#....................................................................................................................................................#
	my ($fastaHsh_ref) = @_;

	my $fastaLengthHsh_ref = {};
	$fastaLengthHsh_ref->{$_} = length ($fastaHsh_ref->{$_}) foreach (keys %{$fastaHsh_ref});

	return ($fastaLengthHsh_ref);
}
sub generateMASTBackgroundFile {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: calculateBackgroundNucleotideFrequency|374, scanMotifWholeGenomeWithMAST|2454
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 5_retrieveSequenceSurroundingPredefinedTSS|213, 6_scanMotifOccurence|235
#	input: $bfilePath, $maxPolymerSize, $seqHsh_ref
#	output: $freqHsh_ref, $proportionHsh_ref
#	toCall: my ($freqHsh_ref, $proportionHsh_ref) = &generateMASTBackgroundFile($seqHsh_ref, $bfilePath, $maxPolymerSize);
#	calledInLine: 408, 2474
#....................................................................................................................................................#
	my ($seqHsh_ref, $bfilePath, $maxPolymerSize) = @_;
	
	my $freqHsh_ref = {};
	my $proportionHsh_ref = {};

	foreach my $polymerSize (1..$maxPolymerSize) {
		my $regexString = '';
		$regexString .= '{A,C,G,T}' foreach (1..$polymerSize);
		$freqHsh_ref->{$polymerSize}{$_} = 0 foreach (glob $regexString);
		foreach my $ID (keys %{$seqHsh_ref}) {
			my $length = length $seqHsh_ref->{$ID};
			foreach my $pos (0..$length-$polymerSize-1) {
				my $polymer = substr $seqHsh_ref->{$ID}, $pos, $polymerSize;
				next if $polymer =~ m/[^ATGC]/i;
				$freqHsh_ref->{$polymerSize}{$polymer}++;
			}
		}
	}
	
	open (FREQ, ">", "$bfilePath");
	foreach my $polymerSize (sort keys %{$freqHsh_ref}) {
		my $sum = sum values %{$freqHsh_ref->{$polymerSize}};
		foreach my $polymer (sort {$freqHsh_ref->{$polymerSize}{$b} <=> $freqHsh_ref->{$polymerSize}{$a}} keys %{$freqHsh_ref->{$polymerSize}}) {
			my $proportion = sprintf "%.10f", $freqHsh_ref->{$polymerSize}{$polymer}/$sum;
			$proportion = 0.0000000001 if $proportion == 0; #--- stupid MAST doesnt like zero!!!
			$proportionHsh_ref->{$polymerSize}{$polymer} = $proportion;
			print FREQ join "", ((join "\t", ($polymer, $proportion)), "\n");
		}
	}
	close FREQ;
	
	return ($freqHsh_ref, $proportionHsh_ref);
}
sub generateShuffleSeq {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|2322
#	appearInSub: >none
#	primaryAppearInSection: 5_retrieveSequenceSurroundingPredefinedTSS|213
#	secondaryAppearInSection: >none
#	input: $resultFastaDir, $resultStorableDir, $seqAroundSiteHsh_ref
#	output: $shuffleSeqAroundSiteHsh_ref, $shuffleSeqAroundSiteInfoHsh_ref
#	toCall: my ($shuffleSeqAroundSiteHsh_ref, $shuffleSeqAroundSiteInfoHsh_ref) = &generateShuffleSeq($seqAroundSiteHsh_ref, $resultFastaDir, $resultStorableDir);
#	calledInLine: 227
#....................................................................................................................................................#
	my ($seqAroundSiteHsh_ref, $resultFastaDir, $resultStorableDir) = @_;
	
	my $shuffleSeqAroundSiteHsh_ref = {};
	my $shuffleSeqAroundSiteInfoHsh_ref = {};
	
	my $shuffleSeqAroundSiteHshPlsPath = "$resultStorableDir/shuffleSeqAroundSiteHsh.pls";
	my $shuffleSeqAroundSiteInfoHshPlsPath = "$resultStorableDir/shuffleSeqAroundSiteInfoHsh.pls";
	
	if (-s $shuffleSeqAroundSiteHshPlsPath and -s $shuffleSeqAroundSiteInfoHshPlsPath) {

		&reportStatus("shuffleSeqAroundSite storable found. Retrieving", 0, "\n");#->2322

		$shuffleSeqAroundSiteHsh_ref = retrieve($shuffleSeqAroundSiteHshPlsPath);
		$shuffleSeqAroundSiteInfoHsh_ref = retrieve($shuffleSeqAroundSiteInfoHshPlsPath);
	
	} else {

		my @itemAry = qw/upStrm dnStrm full/;
		foreach my $siteType (keys %{$seqAroundSiteHsh_ref}) {
			&reportStatus("Shuffling sequences around $siteType", 0, "\n");#->2322
			my %fastaFHHsh = ();
			foreach my $fullOrUpStrmOrDnStrm (@itemAry) {
				my $fastaPath = "$resultFastaDir/shuffle.$siteType.$fullOrUpStrmOrDnStrm.fasta";
				open $fastaFHHsh{$fullOrUpStrmOrDnStrm}, ">", $fastaPath;
				$shuffleSeqAroundSiteInfoHsh_ref->{$siteType}{$fullOrUpStrmOrDnStrm}{'fastaPath'} = $fastaPath;
				$shuffleSeqAroundSiteInfoHsh_ref->{$siteType}{$fullOrUpStrmOrDnStrm}{'totalSeqNum'} = keys %{$seqAroundSiteHsh_ref->{$siteType}};
			}
			foreach my $mRNAID (keys %{$seqAroundSiteHsh_ref->{$siteType}}) {
				my %tmpAryHsh = ();
				@{$tmpAryHsh{'upStrm'}} = split //, $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'};
				@{$tmpAryHsh{'dnStrm'}} = split //, $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'};
				@{$tmpAryHsh{'upStrm'}} = shuffle(@{$tmpAryHsh{'upStrm'}});
				@{$tmpAryHsh{'dnStrm'}} = shuffle(@{$tmpAryHsh{'dnStrm'}});
				@{$tmpAryHsh{'full'}} = (@{$tmpAryHsh{'upStrm'}}, @{$tmpAryHsh{'dnStrm'}});
				foreach my $fullOrUpStrmOrDnStrm (@itemAry) {
					$shuffleSeqAroundSiteHsh_ref->{$siteType}{$fullOrUpStrmOrDnStrm}{$mRNAID} = join '', @{$tmpAryHsh{$fullOrUpStrmOrDnStrm}};
					$shuffleSeqAroundSiteInfoHsh_ref->{$siteType}{$fullOrUpStrmOrDnStrm}{'length'} = length $shuffleSeqAroundSiteHsh_ref->{$siteType}{$mRNAID}{$fullOrUpStrmOrDnStrm} if not $shuffleSeqAroundSiteInfoHsh_ref->{$siteType}{$fullOrUpStrmOrDnStrm}{'length'};
					print {$fastaFHHsh{$fullOrUpStrmOrDnStrm}} ">$mRNAID\n";
					print {$fastaFHHsh{$fullOrUpStrmOrDnStrm}} $shuffleSeqAroundSiteHsh_ref->{$siteType}{$fullOrUpStrmOrDnStrm}{$mRNAID}."\n";
				}
			}
		}

		store($shuffleSeqAroundSiteHsh_ref, $shuffleSeqAroundSiteHshPlsPath);
		store($shuffleSeqAroundSiteInfoHsh_ref, $shuffleSeqAroundSiteInfoHshPlsPath);
		
	}

	return ($shuffleSeqAroundSiteHsh_ref, $shuffleSeqAroundSiteInfoHsh_ref);
}
sub generateThreadHshWithRandomCntg {
#....................................................................................................................................................#
#	subroutineCategory: multithread
#	dependOnSub: >none
#	appearInSub: getCoverageOfItemRngType_multiStrand|852, predictGenomeWidePolyASite|1903
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 9_predictTSS|271
#	input: $cntgAry_ref, $threadToSpawn
#	output: $randCntgInThreadHsh_ref
#	toCall: my ($randCntgInThreadHsh_ref) = &generateThreadHshWithRandomCntg($threadToSpawn, $cntgAry_ref);
#	calledInLine: 866, 1933
#....................................................................................................................................................#

	my ($threadToSpawn, $cntgAry_ref) = @_;

	my @shuffleCntgAry = shuffle(@{$cntgAry_ref});
	my $threadNum = 1;
	my $randCntgInThreadHsh_ref = {};
	foreach my $cntg (@{$cntgAry_ref}) {
		$threadNum = 1 if $threadNum > $threadToSpawn;
		push @{$randCntgInThreadHsh_ref->{$threadNum}}, $cntg;
		$threadNum++;
	}
	
	return ($randCntgInThreadHsh_ref);

}
sub getCoverageOfItemRngType_multiStrand {
#....................................................................................................................................................#
#	subroutineCategory: coverage
#	dependOnSub: generateThreadHshWithRandomCntg|825, reportStatus|2322
#	appearInSub: plotCovAlongCDSBothStrnd|1628
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 10_analyzePAScore|282
#	input: $dirtnAry_ref, $itemByCntgHsh_ref, $itemInfoHsh_ref, $margin, $maxThread, $pileupStorablePathHsh_ref, $rngType
#	output: $covHsh_ref
#	toCall: my ($covHsh_ref) = &getCoverageOfItemRngType_multiStrand($pileupStorablePathHsh_ref, $itemInfoHsh_ref, $itemByCntgHsh_ref, $maxThread, $rngType, $margin, $dirtnAry_ref);
#	calledInLine: 1647
#....................................................................................................................................................#
	my ($pileupStorablePathHsh_ref, $itemInfoHsh_ref, $itemByCntgHsh_ref, $maxThread, $rngType, $margin, $dirtnAry_ref) = @_;
	
	my @cntgAry = (keys %{$pileupStorablePathHsh_ref});
	my $randCntgInThreadHsh_ref = &generateThreadHshWithRandomCntg($maxThread, \@cntgAry);#->825
	my $cntgProc :shared = 0;
	my %threadHsh = ();

	my %strndCovToGetHsh = ();

	$strndCovToGetHsh{'+'}{'s'} = '+';
	$strndCovToGetHsh{'+'}{'a'} = '-';
	$strndCovToGetHsh{'-'}{'s'} = '-';
	$strndCovToGetHsh{'-'}{'a'} = '+';

	foreach my $threadNum (sort {$a <=> $b} keys %{$randCntgInThreadHsh_ref}) {
		my $cntgAry_ref = \@{$randCntgInThreadHsh_ref->{$threadNum}};
		my $cntgNum = @{$randCntgInThreadHsh_ref->{$threadNum}};
		&reportStatus("$cntgNum cntgs spawned to thread $threadNum", 20, "\r");#->2322

		#---spawn a new thread
		($threadHsh{$threadNum}) = threads->new(#---refer to http://www.perlmonks.org/?node_id=966781
		
			sub {
				my ($cntgAry_ref) = @_;
				my $covInThrHsh_ref = {};
				foreach my $cntg (@{$cntgAry_ref}) {
					$cntgProc++;
					my $cntgCovAry_ref = retrieve($pileupStorablePathHsh_ref->{$cntg});
					next if not exists $itemByCntgHsh_ref->{$cntg};
					foreach my $itemID (keys %{$itemByCntgHsh_ref->{$cntg}}) {
						next if not $itemInfoHsh_ref->{$itemID}{$rngType};

						my @rng = sort {$a <=> $b} @{$itemInfoHsh_ref->{$itemID}{$rngType}};
						unshift @rng, ($rng[0]-1-$margin, $rng[0]-1);
						push @rng, ($rng[-1]+1, $rng[-1]+1+$margin);
	
						my $posCovAry_ref = ();
						#my @deBugAry = ($itemID, $itemInfoHsh_ref->{$itemID}{'strnd'});
						for (my $i=0; $i < $#rng; $i += 2) {
							foreach my $j ($rng[$i]-1..$rng[$i+1]-1) {
								my $pos = $j + 1;
								my %tmpCovHsh = ('+'=>0, '-'=>0);
								#push @deBugAry, $j;
								($tmpCovHsh{'+'}, $tmpCovHsh{'-'}) = split /,/, $cntgCovAry_ref->[$pos] if ($cntgCovAry_ref->[$pos]);
								foreach my $dirtn (@{$dirtnAry_ref}) {
									push @{$posCovAry_ref->{$dirtn}}, $tmpCovHsh{$strndCovToGetHsh{$itemInfoHsh_ref->{$itemID}{'strnd'}}{$dirtn}};
								}
							}
						}

						foreach my $dirtn (@{$dirtnAry_ref}) {
							if ($itemInfoHsh_ref->{$itemID}{'strnd'} eq '-') {
								@{$covInThrHsh_ref->{$itemID}{$dirtn}} = reverse @{$posCovAry_ref->{$dirtn}};
							} else {
								@{$covInThrHsh_ref->{$itemID}{$dirtn}} = @{$posCovAry_ref->{$dirtn}};
							}
						}
					}
					
					&reportStatus("Finished counting items on $cntgProc cntg", 20, "\r");#->2322
				}
				return ($covInThrHsh_ref);
			}
			,($cntgAry_ref)
		);
	}
	
	my $covHsh_ref = {};
	
	while (keys %threadHsh) {
		foreach my $threadNum (keys %threadHsh) {
			my $thr = $threadHsh{$threadNum};
			if (not $thr->is_running()) {
				my ($covInThrHsh_ref) = $thr->join;
				foreach my $itemID (keys %{$covInThrHsh_ref}) {
					foreach my $dirtn (keys %{$covInThrHsh_ref->{$itemID}}) {
						$covHsh_ref->{$itemID}{$dirtn} = $covInThrHsh_ref->{$itemID}{$dirtn};
					}
				}
				delete $threadHsh{$threadNum};
				undef $covInThrHsh_ref;
			}
		}
		sleep 1;
	}

	my $numItem = scalar(keys %{$covHsh_ref});
	my $dirtnStr = join " & ", @{$dirtnAry_ref};
	&reportStatus("$rngType coverage of $numItem items in $dirtnStr direction were stored", 10, "\n");#->2322
	
	return ($covHsh_ref);

}
sub getCtgryGeneInfo {
#....................................................................................................................................................#
#	subroutineCategory: general, gff
#	dependOnSub: reportStatus|2322
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|191
#	secondaryAppearInSection: >none
#	input: $ctgryAry_ref, $geneInfoHsh_ref
#	output: $geneCtgryByCntgHsh_ref, $geneCtgryInfoHsh_ref
#	toCall: my ($geneCtgryInfoHsh_ref, $geneCtgryByCntgHsh_ref) = &getCtgryGeneInfo($geneInfoHsh_ref, $ctgryAry_ref);
#	calledInLine: 208
#....................................................................................................................................................#

	my ($geneInfoHsh_ref, $ctgryAry_ref) = @_;
	
	my $geneCtgryInfoHsh_ref = {};
	my $geneCtgryByCntgHsh_ref = {};
	
	my $ctgryStr = join ",", @{$ctgryAry_ref};

	&reportStatus("Filtering GFF on cgtry $ctgryStr", 0, "\n");#->2322
	
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {
		my $ctgry = $geneInfoHsh_ref->{$geneID}{'ctgry'};
		my $cntg = $geneInfoHsh_ref->{$geneID}{'cntg'};
		if (grep /^$ctgry$/, @{$ctgryAry_ref}) {
			%{$geneCtgryInfoHsh_ref->{$geneID}} = %{$geneInfoHsh_ref->{$geneID}};
			$geneCtgryByCntgHsh_ref->{$cntg}{$geneID}++;
		}
	}
	
	my $numGene = keys %{$geneCtgryInfoHsh_ref};
	
	&reportStatus("$numGene gene filtered on cgtry $ctgryStr", 0, "\n");#->2322
	
	return $geneCtgryInfoHsh_ref, $geneCtgryByCntgHsh_ref;
}
sub getMastGenomeBothStrandHit {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: scanMotifWholeGenomeWithMAST|2454
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_scanMotifOccurence|235
#	input: $fastaLengthHsh_ref, $mastHitLog
#	output: $allHitBySeqHsh_ref
#	toCall: my ($allHitBySeqHsh_ref) = &getMastGenomeBothStrandHit($mastHitLog, $fastaLengthHsh_ref);
#	calledInLine: 2489
#....................................................................................................................................................#

	my ($mastHitLog, $fastaLengthHsh_ref) = @_;
	
	my %hitCountHsh = ();
	my $allHitBySeqHsh_ref = {};
	
	open MASTLOG, "<", $mastHitLog;
	while (<MASTLOG>) {
		next if $_ =~ m/^#/;
		chomp;
		my ($sequence_name, $motif, $hit_start, $hit_end, $score, $hit_pValue) = split / +/;
		my $strnd = $1 if $sequence_name =~ m/([\+\-])$/;
		$sequence_name =~ s/[\+\-]$//;
		#print TMPLOG "$sequence_name\t$strnd\t$hit_start\t$hit_pValue\n";
		if ($strnd eq '-') {
			$hit_start = $fastaLengthHsh_ref->{$sequence_name} - $hit_start + 1;
		}
		$allHitBySeqHsh_ref->{$sequence_name}{$strnd}{$hit_start} = $hit_pValue;
	}
	close MASTLOG;
	
	system ("pigz $mastHitLog");
	
	return ($allHitBySeqHsh_ref);
}
sub getMotifBoundCutoff {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: defineMotifPositionBoundsUsingShuffleBackground|644
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 7_defineMotifBoundary|248
#	input: $abvBkgrdFactor, $extraMargin, $minStreakHsh_ref, $motif, $posPctHsh_ref, $siteSeqRng
#	output: $leftBound, $postionFactorHsh_ref, $rightBound
#	toCall: my ($leftBound, $rightBound, $postionFactorHsh_ref) = &getMotifBoundCutoff($posPctHsh_ref, $minStreakHsh_ref, $abvBkgrdFactor, $extraMargin, $siteSeqRng, $motif);
#	calledInLine: 670
#....................................................................................................................................................#
	my ($posPctHsh_ref, $minStreakHsh_ref, $abvBkgrdFactor, $extraMargin, $siteSeqRng, $motif) = @_;
	
	my %regionPosHsh = ();
	my $regionNum = 0;
	my $postvStreak = 0;
	my $negtvStreak = 0;
	my $withinPostvStreak = 'no';
	my %regionSizeHsh = ();
	
	my $queryPeakPct = 0;
	my $queryPeakPos;
	my @posAry = sort {$a <=> $b} keys %{$posPctHsh_ref};
	my $midPtPos = int(@posAry/2)+4;
	foreach my $pos (@posAry) {
		
		#----skip all poses upstream for DSE
		next if ($motif eq 'DSE' and $pos <= $midPtPos);

		if ($posPctHsh_ref->{$pos}{'query'} > $queryPeakPct) {
			$queryPeakPos = $pos;
			$queryPeakPct = $posPctHsh_ref->{$pos}{'query'};
		}
		
		if ($posPctHsh_ref->{$pos}{'query'} > $abvBkgrdFactor*$posPctHsh_ref->{$pos}{'shuffle'}) {
			$postvStreak++;
			$negtvStreak = 0;
			$regionNum++ if $postvStreak == $minStreakHsh_ref->{'positive'} and $withinPostvStreak eq "no";
			$withinPostvStreak = 'yes' if ($postvStreak >= $minStreakHsh_ref->{'positive'});
		} else {
			$negtvStreak++;
			$postvStreak = 0;
			$withinPostvStreak = 'no' if ($negtvStreak >= $minStreakHsh_ref->{'negative'});
		}

		if ($withinPostvStreak eq 'yes') {
			push @{$regionPosHsh{$regionNum}}, $pos;
			$regionSizeHsh{$regionNum}++;
		}
	}
	
	my $leftBound = my $rightBound = 0;
	my $postionFactorHsh_ref = {};
	
	if (exists $regionPosHsh{1} and $motif ne 'CSE') {
		my @regionAry = (sort {$regionSizeHsh{$b} <=> $regionSizeHsh{$a}} keys %regionSizeHsh);
		my $largestRegion = $regionAry[0];
		$leftBound = ${$regionPosHsh{$largestRegion}}[0] - $minStreakHsh_ref->{'positive'} - $extraMargin;
		$rightBound = ${$regionPosHsh{$largestRegion}}[-1] - $minStreakHsh_ref->{'negative'} + $extraMargin;
		my $peakValue = 0;
		for my $pos ($leftBound-1..$rightBound+1) {
			my $pct = $posPctHsh_ref->{$pos}{'query'};
			my $posToTSS = $pos - $siteSeqRng - 1;
			$pct = 0 if $pos < $leftBound or $pos > $rightBound; #---add end zero for better plotting
			$postionFactorHsh_ref->{$posToTSS} = $pct;
			$peakValue = $pct if $pct > $peakValue;
		}
		#---scale positional factor using peak pct value
		$postionFactorHsh_ref->{$_} = $postionFactorHsh_ref->{$_}/$peakValue foreach (keys %{$postionFactorHsh_ref});
	}

	if ($motif eq 'CSE') {
		$leftBound = $rightBound = $queryPeakPos;
		my $posToTSS = $queryPeakPos - $siteSeqRng - 1;
		$postionFactorHsh_ref->{$posToTSS} = 1; 
		$postionFactorHsh_ref->{$posToTSS+1} = 0; #---add end zero for better plotting
		$postionFactorHsh_ref->{$posToTSS-1} = 0; #---add end zero for better plotting
	}

	return ($leftBound, $rightBound, $postionFactorHsh_ref);
}
sub getPAScoreAroundmRNAReferencePoint {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|2322
#	appearInSub: >none
#	primaryAppearInSection: 10_analyzePAScore|282
#	secondaryAppearInSection: >none
#	input: $genomeWidePolyASitePlsPathHsh_ref, $mRNAByCntgHsh_ref, $mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $siteSeqRng
#	output: $PAScoreRefPtBymRNAHsh_ref, $PAScoreRefPtPlotHsh_ref
#	toCall: my ($PAScoreRefPtPlotHsh_ref, $PAScoreRefPtBymRNAHsh_ref) = &getPAScoreAroundmRNAReferencePoint($mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $genomeWidePolyASitePlsPathHsh_ref, $siteSeqRng, $mRNAByCntgHsh_ref);
#	calledInLine: 290
#....................................................................................................................................................#
	my ($mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $genomeWidePolyASitePlsPathHsh_ref, $siteSeqRng, $mRNAByCntgHsh_ref) = @_;
	
	my $PAScoreRefPtPlotHsh_ref = {};
	my $PAScoreRefPtBymRNAHsh_ref = {};
	
	my %tmpStrndDirtnHsh = ();
	$tmpStrndDirtnHsh{'+'}{'s'} = '+';
	$tmpStrndDirtnHsh{'+'}{'a'} = '-';
	$tmpStrndDirtnHsh{'-'}{'s'} = '-';
	$tmpStrndDirtnHsh{'-'}{'a'} = '+';
	
	foreach my $sumOrCount (qw/sum count/) {
		foreach my $siteType (keys %{$mRNARefPtHsh_ref}) {
			foreach my $mode (qw/allPos mostUpStrmPos/) {
				foreach my $rltvPos (-1*$siteSeqRng..$siteSeqRng-1) {
					foreach my $dirtn (qw/a s/) {
						$PAScoreRefPtPlotHsh_ref->{$siteType}{$mode}{$rltvPos}{$dirtn}{$sumOrCount} = 0;
					}
				}
			}
		}
	}

	foreach my $cntg (keys %{$mRNAByCntgHsh_ref}) {
		&reportStatus("Getting TSS scores around site in $cntg", 20, "\r");#->2322

		my $cntgPAScoreAry_ref = retrieve($genomeWidePolyASitePlsPathHsh_ref->{$cntg});
		
		foreach my $mRNAID (keys %{$mRNAByCntgHsh_ref->{$cntg}}) {
			foreach my $siteType (keys %{$mRNARefPtHsh_ref}) {
				my @upStrmPosRng = my @dnStrmPosRng = ();
				if (exists $mRNARefPtHsh_ref->{$siteType}{$mRNAID}) {
					my $geneStrnd = $mRNAInfoHsh_ref->{$mRNAID}{'strnd'};
					my $refPos = $mRNARefPtHsh_ref->{$siteType}{$mRNAID};
					if ($geneStrnd eq '+') {
						@upStrmPosRng = ($refPos-$siteSeqRng-1..$refPos-1);
						@dnStrmPosRng = ($refPos..$refPos+$siteSeqRng);
					} else {
						@upStrmPosRng = reverse ($refPos+1..$refPos+$siteSeqRng+1);
						@dnStrmPosRng = reverse ($refPos-$siteSeqRng..$refPos);
					}
					
					my @searchRngAry = (@upStrmPosRng, @dnStrmPosRng);
					
					my $rltvPos = -1*$siteSeqRng;
					foreach my $pos (@searchRngAry) {
						my $i = $pos - 1;
						if ($cntgPAScoreAry_ref->[$i]) {
							my %tmpScoreHsh = ();
							($tmpScoreHsh{'+'}, $tmpScoreHsh{'-'}) = split ",", $cntgPAScoreAry_ref->[$i];
							foreach my $dirtn (qw/a s/) {
								my $dirtnStrnd = $tmpStrndDirtnHsh{$geneStrnd}{$dirtn};
								$PAScoreRefPtBymRNAHsh_ref->{$mRNAID}{$siteType}{$dirtn}{$rltvPos} = $tmpScoreHsh{$dirtnStrnd} if $tmpScoreHsh{$dirtnStrnd} > 0;
							}
						}
						$rltvPos++;
					}
					foreach my $dirtn (qw/a s/) {
						if (exists $PAScoreRefPtBymRNAHsh_ref->{$mRNAID}{$siteType}{$dirtn}) {
							#---[12/10/2013 13:46] get sorted rltvPos for the $PAScoreRefPtMode eq 'upStrmMost'
							my @rltvPosAry = sort {$a <=> $b} (keys %{$PAScoreRefPtBymRNAHsh_ref->{$mRNAID}{$siteType}{$dirtn}});
							@rltvPosAry = reverse @rltvPosAry if $dirtn eq 'a'; #---[12/10/2013 13:46] reverse for antisense direction, max rltv pos is the upstream
							
							#---[12/10/2013 13:58] go thr the rltv pos from up to down stream
							foreach my $rltvPos (@rltvPosAry) {
								$PAScoreRefPtPlotHsh_ref->{$siteType}{'allPos'}{$rltvPos}{$dirtn}{'sum'} += $PAScoreRefPtBymRNAHsh_ref->{$mRNAID}{$siteType}{$dirtn}{$rltvPos};
								$PAScoreRefPtPlotHsh_ref->{$siteType}{'allPos'}{$rltvPos}{$dirtn}{'count'}++;
							}
							
							#---[12/10/2013 14:05] get only the most upstream position
							my $mostUpStrmPos = $rltvPosAry[0];
							$PAScoreRefPtPlotHsh_ref->{$siteType}{'mostUpStrmPos'}{$mostUpStrmPos}{$dirtn}{'sum'} += $PAScoreRefPtBymRNAHsh_ref->{$mRNAID}{$siteType}{$dirtn}{$mostUpStrmPos};
							$PAScoreRefPtPlotHsh_ref->{$siteType}{'mostUpStrmPos'}{$mostUpStrmPos}{$dirtn}{'count'}++;
						}
					}
				}
			}
		}
	}
	
	return ($PAScoreRefPtPlotHsh_ref, $PAScoreRefPtBymRNAHsh_ref);
}
sub getPAScoreInExonAndTSS {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|2322
#	appearInSub: >none
#	primaryAppearInSection: 10_analyzePAScore|282
#	secondaryAppearInSection: >none
#	input: $genomeWidePolyASitePlsPathHsh_ref, $mRNAByCntgHsh_ref, $mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $siteSeqRng
#	output: $PAScoreExonpolyASiteHshNonZero_ref, $PAScoreExonpolyASiteHshWithZero_ref
#	toCall: my ($PAScoreExonpolyASiteHshNonZero_ref, $PAScoreExonpolyASiteHshWithZero_ref) = &getPAScoreInExonAndTSS($genomeWidePolyASitePlsPathHsh_ref, $mRNAByCntgHsh_ref, $mRNARefPtHsh_ref, $siteSeqRng, $mRNAInfoHsh_ref);
#	calledInLine: 287
#....................................................................................................................................................#
	my ($genomeWidePolyASitePlsPathHsh_ref, $mRNAByCntgHsh_ref, $mRNARefPtHsh_ref, $siteSeqRng, $mRNAInfoHsh_ref) = @_;
	
	
	my $PAScoreExonpolyASiteHshNonZero_ref = {};
	my $PAScoreExonpolyASiteHshWithZero_ref = {};
	
	foreach my $cntg (keys %{$mRNAByCntgHsh_ref}) {
		&reportStatus("Getting polyASite scores at polyASite and exon of genes on $cntg", 20, "\r");#->2322

		my $cntgPAScoreAry_ref = retrieve($genomeWidePolyASitePlsPathHsh_ref->{$cntg});
		
		foreach my $mRNAID (keys %{$mRNAByCntgHsh_ref->{$cntg}}) {
			next if not exists $mRNARefPtHsh_ref->{'polyA_tail'}{$mRNAID};

			my @CDSRngAry = sort {$a <=> $b} @{$mRNAInfoHsh_ref->{$mRNAID}{'CDSRng'}};
			my ($CDSStart, $CDSEnd) = ($CDSRngAry[0], $CDSRngAry[-1]);
			next if $CDSStart+$siteSeqRng > $CDSEnd-$siteSeqRng;

			my %exonPolyASitePosHsh = ();
			@{$exonPolyASitePosHsh{'polyA_tail'}} = ($mRNARefPtHsh_ref->{'polyA_tail'}{$mRNAID});
			
			@{$exonPolyASitePosHsh{'exon'}} = shuffle($CDSStart+$siteSeqRng..$CDSEnd-$siteSeqRng);

			my $geneStrnd = $mRNAInfoHsh_ref->{$mRNAID}{'strnd'};
			foreach my $exonOrPolyASite (keys %exonPolyASitePosHsh) {

				#----will sample non-zero numbers on exon and polyASite, for histogram
				foreach my $pos (@{$exonPolyASitePosHsh{$exonOrPolyASite}}) {
					my $i = $pos-1;
					my %tmpScoreHsh = ();
					$tmpScoreHsh{$exonOrPolyASite}{$_} = 0 foreach (qw/+ -/);
					($tmpScoreHsh{$exonOrPolyASite}{'+'}, $tmpScoreHsh{$exonOrPolyASite}{'-'}) = split ",", $cntgPAScoreAry_ref->[$i] if $cntgPAScoreAry_ref->[$i];
					push @{$PAScoreExonpolyASiteHshNonZero_ref->{$exonOrPolyASite}}, $tmpScoreHsh{$exonOrPolyASite}{$geneStrnd} if $tmpScoreHsh{$exonOrPolyASite}{$geneStrnd} > 0;
				}

				#----will sample a single position on exon and polyASite no matter zero or not for histogram
				my $singlePos = ${$exonPolyASitePosHsh{$exonOrPolyASite}}[0];
				my $i = $singlePos-1;
				my %tmpScoreHsh = ();
				$tmpScoreHsh{$exonOrPolyASite}{$_} = 0 foreach (qw/+ -/);
				($tmpScoreHsh{$exonOrPolyASite}{'+'}, $tmpScoreHsh{$exonOrPolyASite}{'-'}) = split ",", $cntgPAScoreAry_ref->[$i] if $cntgPAScoreAry_ref->[$i];
				push @{$PAScoreExonpolyASiteHshWithZero_ref->{$exonOrPolyASite}}, $tmpScoreHsh{$exonOrPolyASite}{$geneStrnd};
			}
		}
	}

	return ($PAScoreExonpolyASiteHshNonZero_ref, $PAScoreExonpolyASiteHshWithZero_ref);
}
sub getSequenceAroundmRNAReferencePoint {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|2322
#	appearInSub: >none
#	primaryAppearInSection: 5_retrieveSequenceSurroundingPredefinedTSS|213
#	secondaryAppearInSection: >none
#	input: $fastaHsh_ref, $mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $resultFastaDir, $resultStorableDir, $siteSeqRng
#	output: $seqAroundSiteHsh_ref, $seqAroundSiteInfoHsh_ref
#	toCall: my ($seqAroundSiteHsh_ref, $seqAroundSiteInfoHsh_ref) = &getSequenceAroundmRNAReferencePoint($mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $fastaHsh_ref, $siteSeqRng, $resultFastaDir, $resultStorableDir);
#	calledInLine: 222
#....................................................................................................................................................#

	my ($mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $fastaHsh_ref, $siteSeqRng, $resultFastaDir, $resultStorableDir) = @_;
	
	my $seqAroundSiteHsh_ref = {};
	my $seqAroundSiteInfoHsh_ref = {};
	
	my $seqAroundSiteHshPlsPath = "$resultStorableDir/seqAroundSiteHsh.pls";
	my $seqAroundSiteInfoHshPlsPath = "$resultStorableDir/seqAroundSiteInfoHsh.pls";
	
	foreach my $siteType (keys %{$mRNARefPtHsh_ref}) {
		my $fastaPath = "$resultFastaDir/$siteType.full.fasta";

		$seqAroundSiteInfoHsh_ref->{$siteType}{'fastaPath'} = $fastaPath;
		$seqAroundSiteInfoHsh_ref->{$siteType}{'length'} = $siteSeqRng*2;
		$seqAroundSiteInfoHsh_ref->{$siteType}{'totalSeqNum'} = keys %{$mRNARefPtHsh_ref->{$siteType}};
	
		open (FASTA, ">", $fastaPath);
		&reportStatus("Getting sequences around $siteType", 0, "\n");#->2322
		foreach my $mRNAID (keys %{$mRNARefPtHsh_ref->{$siteType}}) {
			my $upStrmSeq;
			my $dnStrmSeq;
			
			if ($mRNAInfoHsh_ref->{$mRNAID}{'strnd'} eq '+') {

				$upStrmSeq = substr $fastaHsh_ref->{$mRNAInfoHsh_ref->{$mRNAID}{'cntg'}}, $mRNARefPtHsh_ref->{$siteType}{$mRNAID}-$siteSeqRng-1, $siteSeqRng;
				$dnStrmSeq= substr $fastaHsh_ref->{$mRNAInfoHsh_ref->{$mRNAID}{'cntg'}}, $mRNARefPtHsh_ref->{$siteType}{$mRNAID}-1, $siteSeqRng;

			} else {
			
				$upStrmSeq = substr $fastaHsh_ref->{$mRNAInfoHsh_ref->{$mRNAID}{'cntg'}}, $mRNARefPtHsh_ref->{$siteType}{$mRNAID}-1+1, $siteSeqRng;
				$dnStrmSeq = substr $fastaHsh_ref->{$mRNAInfoHsh_ref->{$mRNAID}{'cntg'}}, $mRNARefPtHsh_ref->{$siteType}{$mRNAID}-$siteSeqRng-1+1, $siteSeqRng;
			
				$upStrmSeq = reverse $upStrmSeq;
				$upStrmSeq =~ tr/ACGTacgt/TGCAtgca/;

				$dnStrmSeq = reverse $dnStrmSeq;
				$dnStrmSeq =~ tr/ACGTacgt/TGCAtgca/;
			}
			
			if (length($upStrmSeq) == $siteSeqRng and length($dnStrmSeq) == $siteSeqRng) {

				$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'} = $upStrmSeq;
				$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'} = $dnStrmSeq;
			
				print FASTA ">$mRNAID\n";
				print FASTA $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'}.$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'}."\n";
			} else {
				&reportStatus("WARNING:$mRNAID seems to be on cntg edge. Excluding it from sequence analyses", 10, "\n");#->2322
			}
		}
		close FASTA;
	}
	
	store($seqAroundSiteHsh_ref, $seqAroundSiteHshPlsPath);
	store($seqAroundSiteInfoHsh_ref, $seqAroundSiteInfoHshPlsPath);
	
	return ($seqAroundSiteHsh_ref, $seqAroundSiteInfoHsh_ref);
}
sub getSingleMotifMASTLogPostionalData {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: scanMotifAroundSiteWithMAST|2395
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_scanMotifOccurence|235
#	input: $mastHitLog, $maxHitPVal, $maxPos, $totalSeqNum
#	output: $allHitBySeqHsh_ref, $motifPctHsh_ref
#	toCall: my ($motifPctHsh_ref, $allHitBySeqHsh_ref) = &getSingleMotifMASTLogPostionalData($mastHitLog, $maxHitPVal, $maxPos, $totalSeqNum);
#	calledInLine: 2431, 2442
#....................................................................................................................................................#
	my ($mastHitLog, $maxHitPVal, $maxPos, $totalSeqNum) = @_;
	
	my $motifPctHsh_ref = {};
	my %hitCountHsh = ();
	my $allHitBySeqHsh_ref = {};
	
	open MASTLOG, "<", $mastHitLog;
	while (<MASTLOG>) {
		next if $_ =~ m/^#/;
		chomp;
		my ($sequence_name, $motif, $hit_start, $hit_end, $score, $hit_pValue) = split / +/;
		if ($hit_pValue <= $maxHitPVal) {
			$hitCountHsh{$hit_start}++;
		}
		$allHitBySeqHsh_ref->{$sequence_name}{$hit_start} = $hit_pValue;
	}
	close MASTLOG;
	
	foreach my $pos (1..$maxPos) {
		$motifPctHsh_ref->{$pos} = 0;
		$motifPctHsh_ref->{$pos} = 100*$hitCountHsh{$pos}/$totalSeqNum if $hitCountHsh{$pos};
	}

	return ($motifPctHsh_ref, $allHitBySeqHsh_ref);
}
sub getmRNAReferencePoints {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: reportStatus|2322
#	appearInSub: >none
#	primaryAppearInSection: 5_retrieveSequenceSurroundingPredefinedTSS|213
#	secondaryAppearInSection: >none
#	input: $mRNABasedPolyAInfoHsh_ref, $mRNAInfoHsh_ref
#	output: $mRNARefPtHsh_ref
#	toCall: my ($mRNARefPtHsh_ref) = &getmRNAReferencePoints($mRNAInfoHsh_ref, $mRNABasedPolyAInfoHsh_ref);
#	calledInLine: 220
#....................................................................................................................................................#
	my ($mRNAInfoHsh_ref, $mRNABasedPolyAInfoHsh_ref) = @_;

	&reportStatus("Getting mRNA reference point", 10, "\n");#->2322
	my $numGene = 0;
	my $mRNARefPtHsh_ref = {};
	foreach my $mRNAID (sort keys %{$mRNAInfoHsh_ref}) {
		if ($mRNABasedPolyAInfoHsh_ref->{'s'}{'tail'}{$mRNAID}{'peakPos'}) {
			if ($mRNABasedPolyAInfoHsh_ref->{'s'}{'tail'}{$mRNAID}{'peakPos'} ne 'none') {
				$numGene++;
				$mRNARefPtHsh_ref->{"polyA_tail"}{$mRNAID} = $mRNABasedPolyAInfoHsh_ref->{'s'}{'tail'}{$mRNAID}{'peakPos'};
			}
		}
		
		my $strnd = $mRNAInfoHsh_ref->{$mRNAID}{'strnd'};
		my @CDSRngAry = sort {$a <=> $b} @{$mRNAInfoHsh_ref->{$mRNAID}{'CDSRng'}};
		my ($CDSStart, $CDSEnd) = ($CDSRngAry[0], $CDSRngAry[-1]);
		my ($ATGPos, $TAAPos) = ($CDSStart, $CDSEnd);
		if ($strnd eq '-') {
			($ATGPos, $TAAPos) = ($TAAPos, $ATGPos);
			$TAAPos--; #---get the outside position as the reference
		} else {
			$TAAPos++; #---get the outside position as the reference
		}
		
		$mRNARefPtHsh_ref->{'mRNA_ATG'}{$mRNAID} = $ATGPos;
		$mRNARefPtHsh_ref->{'mRNA_TAA'}{$mRNAID} = $TAAPos;
		
	}
	&reportStatus("$numGene mRNA reference point stored", 10, "\n");#->2322
	
	return ($mRNARefPtHsh_ref);
}
sub ggplotMultiSampleBoxWhisker {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: plotPAScoreInExonAndTSS|1805
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 10_analyzePAScore|282
#	input: $RScriptPath, $dataPath, $dataPtMax, $extraArg, $height, $log2OrLinear, $logPath, $pdfPath, $plotAryHsh_ref, $width, $yAxis
#	output: 
#	toCall: &ggplotMultiSampleBoxWhisker($plotAryHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $yAxis, $extraArg, $log2OrLinear, $dataPtMax, $height, $width);
#	calledInLine: 1857
#....................................................................................................................................................#
	my ($plotAryHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $yAxis, $extraArg, $log2OrLinear, $dataPtMax, $height, $width) = @_;
	
	open (PLOTDATA, ">", $dataPath);
	print PLOTDATA "sample\t$yAxis\n";

	foreach my $sample (keys %{$plotAryHsh_ref}) {
		my $plotAry_ref = $plotAryHsh_ref->{$sample};

		my $indivDataPtMax = $dataPtMax;
		$indivDataPtMax = @{$plotAry_ref} if $indivDataPtMax > @{$plotAry_ref};

		#---down sample the data point number
		my @shuffleIndexAry = shuffle(0..$#{$plotAry_ref});
		foreach my $i (0..$indivDataPtMax-1) {
			my $value = $plotAry_ref->[$shuffleIndexAry[$i]];
			print PLOTDATA "$sample\t$value\n";
		}
		print PLOTDATA "\n";
	}
	
	my $scale = '';
	$scale = ' + scale_y_continuous(trans=log2_trans())' if $log2OrLinear eq 'log2';
	open (R, ">", $RScriptPath);
	print R "library(ggplot2)"."\n";
	print R "library(scales)"."\n";
	print R "dataFrame = read.table(file=\"$dataPath\", header=TRUE)"."\n";
	print R "ggplot(dataFrame, aes(x=sample, y=$yAxis)) + ggtitle(\"Boxplot of $yAxis in $log2OrLinear scale\") + geom_boxplot(aes(fill = sample)) $scale $extraArg"."\n";
	print R "ggsave(file=\"$pdfPath\", height=$height, width=$width)\n";
	
	system ("R --slave --vanilla --file=$RScriptPath &>$logPath ");

	return ();
}
sub ggplotTwoSampleHistogram {
#....................................................................................................................................................#
#	subroutineCategory: ggplot, plotInR
#	dependOnSub: >none
#	appearInSub: plotPAScoreInExonAndTSS|1805
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 10_analyzePAScore|282
#	input: $RScriptPath, $binWidth, $dataPath, $dataPtMax, $densityOrFrequency, $extraArg, $leftxAxisPercentileLimit, $log2OrLinear, $logPath, $pdfPath, $plotAryHsh_ref, $rightxAxisPercentileLimit, $xAxis
#	output: none
#	toCall: &ggplotTwoSampleHistogram($plotAryHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $leftxAxisPercentileLimit, $rightxAxisPercentileLimit, $xAxis, $binWidth, $dataPtMax, $extraArg, $log2OrLinear, $densityOrFrequency);
#	calledInLine: 1843
#....................................................................................................................................................#

	my ($plotAryHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $leftxAxisPercentileLimit, $rightxAxisPercentileLimit, $xAxis, $binWidth, $dataPtMax, $extraArg, $log2OrLinear, $densityOrFrequency) = @_;
	
	open (PLOTDATA, ">", $dataPath);
	print PLOTDATA "sample\t$xAxis\n";

	foreach my $sample (keys %{$plotAryHsh_ref}) {
		
		my $plotAry_ref = $plotAryHsh_ref->{$sample};
	
		my $valueStatObj = Statistics::Descriptive::Full->new();
		$valueStatObj->add_data(@{$plotAry_ref});

		my $leftxAxisLimitValue;
		if ($leftxAxisPercentileLimit eq 'min') {
			$leftxAxisLimitValue = $valueStatObj->min();
		} else {
			$leftxAxisLimitValue = $valueStatObj->percentile($leftxAxisPercentileLimit);
		}

		my $rightxAxisLimitValue;
		if ($rightxAxisPercentileLimit eq 'max') {
			$rightxAxisLimitValue = $valueStatObj->max();
		} else {
			$rightxAxisLimitValue = $valueStatObj->percentile($rightxAxisPercentileLimit);
		}
	
		#---trim the end values
		my @trimmedAry = ();
		foreach my $value (@{$plotAry_ref}) {
			push @trimmedAry, $value if $value <= $rightxAxisLimitValue and $value >= $leftxAxisLimitValue;
		}
		
		my $indivDataPtMax = $dataPtMax;
		$indivDataPtMax = @trimmedAry if $indivDataPtMax > @trimmedAry;

		#---down sample the data point number
		my @shuffleIndexAry = shuffle(0..$#trimmedAry);
		foreach my $i (0..$indivDataPtMax-1) {
			my $shuffleValue = $trimmedAry[$shuffleIndexAry[$i]];
			print PLOTDATA "$sample\t$shuffleValue\n";
		}
	}
	close PLOTDATA;

	my $scale = '';
	$scale = ' + scale_x_continuous(trans=log2_trans())' if $log2OrLinear eq 'log2';

	my $plotGeom = "ggplot(dataFrame, aes($xAxis, fill = sample)) + geom_density(alpha = 0.2)";#---default = density
	$plotGeom = "ggplot(dataFrame, aes($xAxis, color = sample)) + geom_freqpoly(binwidth = $binWidth)" if $densityOrFrequency eq 'frequency';

	open (R, ">", $RScriptPath);
	print R "library(ggplot2)"."\n";
	print R "library(scales)"."\n";
	print R "dataFrame = read.table(file=\"$dataPath\", header=TRUE)"."\n";
	print R "$plotGeom + ggtitle(\"Distribution of $xAxis $log2OrLinear scale\") $scale $extraArg"."\n";
	print R "ggsave(file=\"$pdfPath\")\n";
	close R;
	
	system ("R --slave --vanilla --file=$RScriptPath &>$logPath ");

}
sub ggplotXYLinesMultipleSamples {
#....................................................................................................................................................#
#	subroutineCategory: ggplot, plotInR
#	dependOnSub: >none
#	appearInSub: defineMotifPositionBoundsUsingShuffleBackground|644, plotBaseCompositionAroundmRNAReferencePoint|1572, plotCovAlongCDSBothStrnd|1628, plotMotifPostionFactor|1724, plotPAScoreAroundmRNAReferencePoint|1761
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 10_analyzePAScore|282, 5_retrieveSequenceSurroundingPredefinedTSS|213, 7_defineMotifBoundary|248
#	input: $RScriptPath, $XAXis, $YAxis, $YVariable, $dataPath, $extraArg, $height, $logPath, $pdfPath, $plotDataHsh_ref, $width
#	output: none
#	toCall: &ggplotXYLinesMultipleSamples($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $XAXis, $YAxis, $YVariable, $extraArg, $height, $width);
#	calledInLine: 692, 1623, 1718, 1756, 1797
#....................................................................................................................................................#

	my ($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $XAXis, $YAxis, $YVariable, $extraArg, $height, $width) = @_;
	
	open (PLOTDATA, ">", $dataPath);
	print PLOTDATA join "", (join "\t", ($YVariable, $YAxis, $XAXis)), "\n";
	foreach my $YCategory (sort keys %{$plotDataHsh_ref}) {
		foreach my $XVal (sort {$a <=> $b} keys %{$plotDataHsh_ref->{$YCategory}}) {
			my $YVal = $plotDataHsh_ref->{$YCategory}{$XVal};
			print PLOTDATA join "", (join "\t", ($YCategory, $YVal, $XVal)), "\n";
		}
	}
	close PLOTDATA;

	open (R, ">", $RScriptPath);
	print R "library(ggplot2)"."\n";
	print R "dataFrame = read.table(file=\"$dataPath\", header=TRUE)"."\n";
	print R "ggplot(dataFrame, aes(x=$XAXis, y=$YAxis, colour=$YVariable)) + geom_line() $extraArg"."\n";
	print R "ggsave(file=\"$pdfPath\", height=$height, width=$width)\n";
	close R;
	
	system ("R --slave --vanilla --file=$RScriptPath &>$logPath");

}
sub plotBaseCompositionAroundmRNAReferencePoint {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: calculateBaseCompositionInAlignments|420, ggplotXYLinesMultipleSamples|1537, reportStatus|2322
#	appearInSub: >none
#	primaryAppearInSection: 5_retrieveSequenceSurroundingPredefinedTSS|213
#	secondaryAppearInSection: >none
#	input: $ggplotDirHsh_ref, $seqAroundSiteHsh_ref
#	output: none
#	toCall: &plotBaseCompositionAroundmRNAReferencePoint($seqAroundSiteHsh_ref, $ggplotDirHsh_ref);
#	calledInLine: 224
#....................................................................................................................................................#

	my ($seqAroundSiteHsh_ref, $ggplotDirHsh_ref) = @_;
	
	my $tmpItemHsh_ref = {};
	$tmpItemHsh_ref->{'polyA_tail'} = 'baseComp_polyA_tail';
	$tmpItemHsh_ref->{'mRNA_ATG'} = 'baseComp_mRNA_ATG';
	$tmpItemHsh_ref->{'mRNA_TAA'} = 'baseComp_mRNA_TAA';

	my $sub_ggplotDir = "baseCompstn";

	system "mkdir -pm 777 $ggplotDirHsh_ref->{$_}/$sub_ggplotDir/" foreach (keys %{$ggplotDirHsh_ref});
	
	foreach my $siteType (sort keys %{$seqAroundSiteHsh_ref}) {
		
		&reportStatus("Plotting sequence compostion around $siteType", 10, "\n");#->2322

		my $item = $tmpItemHsh_ref->{$siteType};
		
		my $seqAlignHsh_ref = {};
		foreach my $mRNAID (sort keys %{$seqAroundSiteHsh_ref->{$siteType}}) {
			my $fullSeq = $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'}.$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'};
			$seqAlignHsh_ref->{$mRNAID} = $fullSeq;
		}


		my $seqLength = length ($seqAlignHsh_ref->{(keys %{$seqAlignHsh_ref})[rand keys %{$seqAlignHsh_ref}]}); #---from my $random_value = $hash{(keys %hash)[rand keys %hash]} in http://stackoverflow.com/questions/8547642/select-a-random-hash-key
		my $baseCompByBaseHsh_ref = &calculateBaseCompositionInAlignments($seqAlignHsh_ref);#->420
		
		{	
			my $plotDataHsh_ref = $baseCompByBaseHsh_ref;
			my $dataPath = "$ggplotDirHsh_ref->{'dat'}/$sub_ggplotDir/$item.dat";
			my $pdfPath = "$ggplotDirHsh_ref->{'pdf'}/$sub_ggplotDir/$item.pdf";
			my $RScriptPath = "$ggplotDirHsh_ref->{'R'}/$sub_ggplotDir/$item.R";
			my $logPath = "$ggplotDirHsh_ref->{'log'}/$sub_ggplotDir/$item.log";
			my $XAXis = 'relativePositon';
			my $YAxis = 'proportion';
			my $YVariable = 'base';
			my $extraArg = "+ ylim(0,1) + scale_x_continuous(breaks=seq(0, $seqLength, by=5))";
			my $height = 9;
			my $width = 18;
			&ggplotXYLinesMultipleSamples($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $XAXis, $YAxis, $YVariable, $extraArg, $height, $width);#->1537
		}
	}
}
sub plotCovAlongCDSBothStrnd {
#....................................................................................................................................................#
#	subroutineCategory: coverage
#	dependOnSub: getCoverageOfItemRngType_multiStrand|852, ggplotXYLinesMultipleSamples|1537, reportStatus|2322
#	appearInSub: >none
#	primaryAppearInSection: 10_analyzePAScore|282
#	secondaryAppearInSection: >none
#	input: $binWidth, $ggplotDirHsh_ref, $mRNAByCntgHsh_ref, $mRNAInfoHsh_ref, $maxThread, $minCDSLength, $pileupStorablePathHsh_ref
#	output: $covAlongCDSBymRNAHsh_ref, $pooledCovAlongCDSHsh_ref
#	toCall: my ($pooledCovAlongCDSHsh_ref, $covAlongCDSBymRNAHsh_ref) = &plotCovAlongCDSBothStrnd($pileupStorablePathHsh_ref, $mRNAInfoHsh_ref, $mRNAByCntgHsh_ref, $minCDSLength, $ggplotDirHsh_ref, $maxThread, $binWidth);
#	calledInLine: 294
#....................................................................................................................................................#
	my ($pileupStorablePathHsh_ref, $mRNAInfoHsh_ref, $mRNAByCntgHsh_ref, $minCDSLength, $ggplotDirHsh_ref, $maxThread, $binWidth) = @_;
	
	#---get the coverage
	my $itemInfoHsh_ref = $mRNAInfoHsh_ref;
	my $itemByCntgHsh_ref = $mRNAByCntgHsh_ref;
	my $rngType = 'CDSRng';
	my $margin = -1;
	my $dirtnAry_ref = [qw/a s/];
	my ($covHsh_ref) = &getCoverageOfItemRngType_multiStrand($pileupStorablePathHsh_ref, $itemInfoHsh_ref, $itemByCntgHsh_ref, $maxThread, $rngType, $margin, $dirtnAry_ref);#->852
	
	my $pooledCovAlongCDSHsh_ref = {};
	my $covAlongCDSBymRNAHsh_ref = {};
	
	foreach my $ATGOrTAA (qw/ATG TAA/) {
		foreach my $dirtn (qw/a s/) {
			foreach my $rltvPos (1..$minCDSLength) {
				$pooledCovAlongCDSHsh_ref->{$ATGOrTAA}{$dirtn}{$rltvPos} = 0;
			}
		}
	}

	&reportStatus("Collecting mRNA coverage along CDS", 10, "\n");#->2322
	
	my $numGene = 0;
	
	foreach my $mRNAID (sort keys %{$mRNAInfoHsh_ref}) {
		my @CDSRngAry = sort {$a <=> $b} @{$mRNAInfoHsh_ref->{$mRNAID}{'CDSRng'}};
		my $CDSLength = $CDSRngAry[-1] - $CDSRngAry[0] + 1;
		
		#---[15/10/2013 18:05] only genes with longer than minCDSLength will be plot 
		if ($CDSLength >= $minCDSLength) {
			$numGene++;
			foreach my $ATGOrTAA (qw/ATG TAA/) {
				foreach my $dirtn (qw/a s/) {
					my @covAry = ();
					if ($ATGOrTAA eq 'ATG') {
						#---[15/10/2013 18:01] start from ATG
						@covAry = @{$covHsh_ref->{$mRNAID}{$dirtn}};
					} elsif ($ATGOrTAA eq 'TAA') {
						#---[15/10/2013 18:01] start from TAA
						@covAry = reverse @{$covHsh_ref->{$mRNAID}{$dirtn}};
					}
					my $rltvPos = 0;
					
					my $upStrmMostRecorded = 'no';
					
					foreach my $cov (@covAry) {
						$rltvPos++;
						$cov = 0 if $upStrmMostRecorded eq 'yes';
						$pooledCovAlongCDSHsh_ref->{$ATGOrTAA}{$dirtn}{$rltvPos} += $cov;
						$covAlongCDSBymRNAHsh_ref->{$ATGOrTAA}{$dirtn}{$mRNAID}{$rltvPos} = $cov;
						$upStrmMostRecorded = 'yes' if $cov > 0;
						last if $rltvPos >= $minCDSLength;
					}
				}
			}
		}
	}
	
	foreach my $ATGOrTAA (keys %{$pooledCovAlongCDSHsh_ref}) {
		my $plotDataHsh_ref = {};
		foreach my $dirtn (keys %{$pooledCovAlongCDSHsh_ref->{$ATGOrTAA}}) {
			foreach my $rltvPos (keys %{$pooledCovAlongCDSHsh_ref->{$ATGOrTAA}{$dirtn}}) {
				my $binDistance  = nearest($binWidth, $rltvPos);
				$plotDataHsh_ref->{$dirtn}{$binDistance} = 0 if not $plotDataHsh_ref->{$dirtn}{$binDistance};
				$plotDataHsh_ref->{$dirtn}{$binDistance} += $pooledCovAlongCDSHsh_ref->{$ATGOrTAA}{$dirtn}{$rltvPos};
			}
		}
		my $nameTag = "PAScore.from.$ATGOrTAA.alongCDS";
		my $pdfPath = $ggplotDirHsh_ref->{'pdf'}."/$nameTag.pdf";
		my $dataPath = $ggplotDirHsh_ref->{'dat'}."/$nameTag.dat";
		my $RScriptPath = $ggplotDirHsh_ref->{'R'}."/$nameTag.R";
		my $logPath = $ggplotDirHsh_ref->{'log'}."/$nameTag.log";
		my $xAxis = 'relative_positon';
		my $YAxis = 'PAScore';
		my $YVariable = 'direction';
		my $extraArg = "+ ggtitle(\"N=$numGene minCDSLength=$minCDSLength\")";
		my $height = 6;
		my $width = 14;
		&ggplotXYLinesMultipleSamples($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $xAxis, $YAxis, $YVariable, $extraArg, $height, $width);#->1537
	}
	
	return ($pooledCovAlongCDSHsh_ref, $covAlongCDSBymRNAHsh_ref);
}
sub plotMotifPostionFactor {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: ggplotXYLinesMultipleSamples|1537
#	appearInSub: >none
#	primaryAppearInSection: 7_defineMotifBoundary|248
#	secondaryAppearInSection: >none
#	input: $generalggplotDirHsh_ref, $motifPostionFactorHsh_ref
#	output: 
#	toCall: &plotMotifPostionFactor($motifPostionFactorHsh_ref, $generalggplotDirHsh_ref);
#	calledInLine: 254
#....................................................................................................................................................#
	my ($motifPostionFactorHsh_ref, $generalggplotDirHsh_ref) = @_;
	
	my $plotDataHsh_ref = {};
	foreach my $motif (keys %{$motifPostionFactorHsh_ref}) {
		foreach my $pos (keys %{$motifPostionFactorHsh_ref->{$motif}}) {
			$plotDataHsh_ref->{$motif}{$pos} = $motifPostionFactorHsh_ref->{$motif}{$pos};
		}
	}
	
	my $nameTag = "motifPostionFactor";
	my $pdfPath = $generalggplotDirHsh_ref->{'pdf'}."/$nameTag.pdf";
	my $dataPath = $generalggplotDirHsh_ref->{'dat'}."/$nameTag.dat";
	my $RScriptPath = $generalggplotDirHsh_ref->{'R'}."/$nameTag.R";
	my $logPath = $generalggplotDirHsh_ref->{'log'}."/$nameTag.log";
	my $xAxis = 'relativePositon';
	my $YAxis = 'factor';
	my $YVariable = 'motif';
	#my $extraArg = '+ ylim(0, 100)';
	my $extraArg = '';
	my $height = 6;
	my $width = 14;
	&ggplotXYLinesMultipleSamples($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $xAxis, $YAxis, $YVariable, $extraArg, $height, $width);#->1537

	return ();
}
sub plotPAScoreAroundmRNAReferencePoint {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: ggplotXYLinesMultipleSamples|1537
#	appearInSub: >none
#	primaryAppearInSection: 10_analyzePAScore|282
#	secondaryAppearInSection: >none
#	input: $PAScoreRefPtPlotHsh_ref, $PAScoreggplotDirHsh_ref, $mRNARefPtHsh_ref
#	output: 
#	toCall: &plotPAScoreAroundmRNAReferencePoint($PAScoreRefPtPlotHsh_ref, $PAScoreggplotDirHsh_ref, $mRNARefPtHsh_ref);
#	calledInLine: 292
#....................................................................................................................................................#
	my ($PAScoreRefPtPlotHsh_ref, $PAScoreggplotDirHsh_ref, $mRNARefPtHsh_ref) = @_;
	
	#---[12/10/2013 16:42] choose to plot the s=um of the TSSSocre or just the frequency count
	#foreach my $sumOrCount (qw/sum count/) {
	foreach my $sumOrCount (qw/sum/) {
		foreach my $siteType (keys %{$PAScoreRefPtPlotHsh_ref}) {
			my $totalNum = keys %{$mRNARefPtHsh_ref->{$siteType}};
			foreach my $mode (keys %{$PAScoreRefPtPlotHsh_ref->{$siteType}}) {
				my $plotDataHsh_ref = {};
				foreach my $rltvPos (sort keys %{$PAScoreRefPtPlotHsh_ref->{$siteType}{$mode}}) {
					foreach my $dirtn (sort keys %{$PAScoreRefPtPlotHsh_ref->{$siteType}{$mode}{$rltvPos}}) {
						$plotDataHsh_ref->{$dirtn}{$rltvPos} = $PAScoreRefPtPlotHsh_ref->{$siteType}{$mode}{$rltvPos}{$dirtn}{$sumOrCount};
					}
				}
				my $nameTag = "$sumOrCount.$mode.$siteType.TSS.Score";
				my $pdfPath = $PAScoreggplotDirHsh_ref->{'pdf'}."/$nameTag.pdf";
				my $dataPath = $PAScoreggplotDirHsh_ref->{'dat'}."/$nameTag.dat";
				my $RScriptPath = $PAScoreggplotDirHsh_ref->{'R'}."/$nameTag.R";
				my $logPath = $PAScoreggplotDirHsh_ref->{'log'}."/$nameTag.log";
				my $xAxis = 'relative_Positon';
				my $YAxis = 'average_TSS_Score';
				my $YVariable = 'direction';
				my $extraArg = "+ ggtitle(\"N=$totalNum\")";
				my $height = 6;
				my $width = 14;
				&ggplotXYLinesMultipleSamples($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $xAxis, $YAxis, $YVariable, $extraArg, $height, $width);#->1537
			}
		}
	}
	
	return ();
}
sub plotPAScoreInExonAndTSS {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: ggplotMultiSampleBoxWhisker|1418, ggplotTwoSampleHistogram|1463
#	appearInSub: >none
#	primaryAppearInSection: 10_analyzePAScore|282
#	secondaryAppearInSection: >none
#	input: $PAScoreExonTSSHshNonZero_ref, $PAScoreExonTSSHshWithZero_ref, $PAScoreggplotDirHsh_ref
#	output: 
#	toCall: &plotPAScoreInExonAndTSS($PAScoreExonTSSHshNonZero_ref, $PAScoreExonTSSHshWithZero_ref, $PAScoreggplotDirHsh_ref);
#	calledInLine: 288
#....................................................................................................................................................#
	my ($PAScoreExonTSSHshNonZero_ref, $PAScoreExonTSSHshWithZero_ref, $PAScoreggplotDirHsh_ref) = @_;
	
	my %tmpHsh = ();
	$tmpHsh{'withZero'}{'plotAryHsh_ref'} = $PAScoreExonTSSHshWithZero_ref;
	$tmpHsh{'nonZero'}{'plotAryHsh_ref'} = $PAScoreExonTSSHshNonZero_ref;
	$tmpHsh{'withZero'}{'densityOrFrequency'} = 'density';
	$tmpHsh{'nonZero'}{'densityOrFrequency'} = 'density';
	$tmpHsh{'withZero'}{'log2OrLinear'} = 'linear';
	$tmpHsh{'nonZero'}{'log2OrLinear'} = 'linear';
	
	foreach my $withZeroOrNonZero (keys %tmpHsh) {
		my $plotAryHsh_ref = $tmpHsh{$withZeroOrNonZero}{'plotAryHsh_ref'};
		my $dataPtMax = 3000;

		{
			my $nameTag = "PAScore.TSS.vs.Exon.histogram.$withZeroOrNonZero";
			my $dataPath = "$PAScoreggplotDirHsh_ref->{dat}/$nameTag.dat";
			my $pdfPath = "$PAScoreggplotDirHsh_ref->{pdf}/$nameTag.pdf";
			my $RScriptPath = "$PAScoreggplotDirHsh_ref->{R}/$nameTag.R";
			my $logPath = "$PAScoreggplotDirHsh_ref->{log}/$nameTag.log";
			my $leftxAxisPercentileLimit = 'min';
			my $rightxAxisPercentileLimit = 'max';
			my $xAxis = "PAScore";
			my $binWidth = 0.1;
			my $log2OrLinear = $tmpHsh{$withZeroOrNonZero}{'log2OrLinear'};
			my $extraArg = '';
			my $densityOrFrequency = $tmpHsh{$withZeroOrNonZero}{'densityOrFrequency'};
			&ggplotTwoSampleHistogram($plotAryHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $leftxAxisPercentileLimit, $rightxAxisPercentileLimit, $xAxis, $binWidth, $dataPtMax, $extraArg, $log2OrLinear, $densityOrFrequency);#->1463
		}
	
		{
			my $nameTag = "PAScore.TSS.vs.Exon.box.$withZeroOrNonZero";
			my $dataPath = "$PAScoreggplotDirHsh_ref->{dat}/$nameTag.dat";
			my $pdfPath = "$PAScoreggplotDirHsh_ref->{pdf}/$nameTag.pdf";
			my $RScriptPath = "$PAScoreggplotDirHsh_ref->{R}/$nameTag.R";
			my $logPath = "$PAScoreggplotDirHsh_ref->{log}/$nameTag.log";
			my $yAxis = "PAScore";
			my $log2OrLinear = $tmpHsh{$withZeroOrNonZero}{'log2OrLinear'};
			my $extraArg = '';
			my $height = 12;
			my $width = 8;
			&ggplotMultiSampleBoxWhisker($plotAryHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $yAxis, $extraArg, $log2OrLinear, $dataPtMax, $height, $width);#->1418
		}
	}
	return ();
}
sub plotWeblogoAroundRefPoint {
#....................................................................................................................................................#
#	subroutineCategory: specific, thirdPartyApp
#	dependOnSub: createWeblogo|600
#	appearInSub: >none
#	primaryAppearInSection: 5_retrieveSequenceSurroundingPredefinedTSS|213
#	secondaryAppearInSection: >none
#	input: $seqAroundSiteHsh_ref, $weblogoDirHsh_ref
#	output: none
#	toCall: &plotWeblogoAroundRefPoint($seqAroundSiteHsh_ref, $weblogoDirHsh_ref);
#	calledInLine: 225
#....................................................................................................................................................#
	
	my ($seqAroundSiteHsh_ref, $weblogoDirHsh_ref) = @_;
	
	my $upStrmRng = 10;
	my $dnStrmRng = 10;
	
	foreach my $siteType (keys %{$seqAroundSiteHsh_ref}) {
		my $seqAlignHsh_ref = {};
		foreach my $mRNAID (keys %{$seqAroundSiteHsh_ref->{$siteType}}) {
			my $upStrmSeq = substr $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'}, -1*$upStrmRng, $upStrmRng;
			my $dnStrmRng = substr $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'}, 0, $dnStrmRng;
			$seqAlignHsh_ref->{'sense'}{$mRNAID} = $upStrmSeq.$dnStrmRng;
			$seqAlignHsh_ref->{'antisense'}{$mRNAID} = reverse $seqAlignHsh_ref->{'sense'}{$mRNAID};
			$seqAlignHsh_ref->{'antisense'}{$mRNAID} =~ tr/ACGTacgt/TGCAtgca/;
		}
		
		foreach my $dirtn (keys %{$seqAlignHsh_ref}) {
			my $numSeq = keys %{$seqAroundSiteHsh_ref->{$siteType}};
			my $nameTag = "around.$siteType.$dirtn";
			my $pdfPath = "$weblogoDirHsh_ref->{pdf}/$nameTag.$dirtn.pdf";
			my $fastaPath = "$weblogoDirHsh_ref->{fasta}/$nameTag.$dirtn.fasta";
			my $title = "$siteType\_-$upStrmRng..+$dnStrmRng\[n=$numSeq\]";
			my $seqType = 'dna';
			&createWeblogo($seqAlignHsh_ref->{$dirtn}, $pdfPath, $fastaPath, $seqType, $title);#->600
		}
	}
	
}
sub predictGenomeWidePolyASite {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: checkRunningThreadAndWaitToJoin|490, createEmptyStorableForGenowideTSSPredictionData|518, generateThreadHshWithRandomCntg|825, reportStatus|2322
#	appearInSub: >none
#	primaryAppearInSection: 9_predictTSS|271
#	secondaryAppearInSection: >none
#	input: $PAScoreStorableDir, $cntgMotifHitPlsPathHsh_ref, $fastaHsh_ref, $logTransformPVal, $motifPostionFactorHsh_ref, $predictMotifInfoHsh_ref
#	output: $genomeWidePolyASitePlsPathHsh_ref
#	toCall: my ($genomeWidePolyASitePlsPathHsh_ref) = &predictGenomeWidePolyASite($motifPostionFactorHsh_ref, $cntgMotifHitPlsPathHsh_ref, $predictMotifInfoHsh_ref, $fastaHsh_ref, $PAScoreStorableDir, $logTransformPVal);
#	calledInLine: 276
#....................................................................................................................................................#
	my ($motifPostionFactorHsh_ref, $cntgMotifHitPlsPathHsh_ref, $predictMotifInfoHsh_ref, $fastaHsh_ref, $PAScoreStorableDir, $logTransformPVal) = @_;
	
	my ($genomeWidePolyASitePlsPathHsh_ref, $allGenomeWideTSSPlsPathExist) = &createEmptyStorableForGenowideTSSPredictionData($fastaHsh_ref, $PAScoreStorableDir);#->518
	
	#----calculate the TSS Score only if not all GenomeWideTSSPlsPath exists
	if ($allGenomeWideTSSPlsPathExist eq 'no') {
		&reportStatus("Start calculating the TSS Score in all cntg", 0, "\n");#->2322

		my $motifInfoHsh_ref = {};
		foreach my $motif (keys %{$motifPostionFactorHsh_ref}) {
			$motifInfoHsh_ref->{$motif}{'postionFactorHsh_ref'} = $motifPostionFactorHsh_ref->{$motif};
			@{$motifInfoHsh_ref->{$motif}{'rltvPosAry'}} = (sort {$a <=> $b} keys %{$motifPostionFactorHsh_ref->{$motif}});
			$motifInfoHsh_ref->{$motif}{'mustValid'} = $predictMotifInfoHsh_ref->{$motif}{'mustValid'};
			$motifInfoHsh_ref->{$motif}{'maxPValGenomePredict'} = $predictMotifInfoHsh_ref->{$motif}{'maxPValGenomePredict'};
			$motifInfoHsh_ref->{$motif}{'scoreFactor'} = $predictMotifInfoHsh_ref->{$motif}{'scoreFactor'};
		}
	
		my $threadToSpawn = 10;
		my @cntgAry = (keys %{$genomeWidePolyASitePlsPathHsh_ref});
		my $randCntgInThreadHsh_ref = &generateThreadHshWithRandomCntg($threadToSpawn, \@cntgAry);#->825
		my $cntgProc :shared = 0;
		foreach my $threadNum (sort {$a <=> $b} keys %{$randCntgInThreadHsh_ref}) {
			my $cntgAry_ref = \@{$randCntgInThreadHsh_ref->{$threadNum}};
			my $cntgNum = @{$randCntgInThreadHsh_ref->{$threadNum}};

			&reportStatus("$cntgNum cntgs spawned to thread $threadNum", 0, "\n");#->2322

			#---spawn a new thread
			threads->create(
		
				sub {
					my ($cntgAry_ref) = @_;

						foreach my $cntg (@{$cntgAry_ref}) {
							$cntgProc++;
	
							&reportStatus("$cntgProc cntgs processed", 20, "\r");#->2322

							my $cntgPAScoreAry_ref = retrieve($genomeWidePolyASitePlsPathHsh_ref->{$cntg});
							my $cntgMotifHitHsh_ref = retrieve($cntgMotifHitPlsPathHsh_ref->{$cntg});
						
							foreach my $i (0..$#{$cntgPAScoreAry_ref}) {
								my $pos = $i + 1;
								my %strndScoreHsh = ();
								foreach my $strnd (qw/+ -/) {
									$strndScoreHsh{$strnd} = 0;
									foreach my $motif (keys %{$motifInfoHsh_ref}) {
										my %motifScoreHsh = ();
										$motifScoreHsh{$motif} = 0;
										my $postionFactorHsh_ref = $motifInfoHsh_ref->{$motif}{'postionFactorHsh_ref'};
					
										foreach my $rltvPos (@{$motifInfoHsh_ref->{$motif}{'rltvPosAry'}}) {
											my $srchPos;
											if ($strnd eq '+') {
												$srchPos = $pos + $rltvPos;
											} else {
												$srchPos = $pos - $rltvPos;
											}
						
											if ($cntgMotifHitHsh_ref->{$strnd}{$srchPos}{$motif}) {
												my $pval = $cntgMotifHitHsh_ref->{$strnd}{$srchPos}{$motif};
												if ($pval <= $motifInfoHsh_ref->{$motif}{'maxPValGenomePredict'}) {
													my $postionFactor = $postionFactorHsh_ref->{$rltvPos};
													my $scoreFactor = $motifInfoHsh_ref->{$motif}{'scoreFactor'};
													my $transformPVal = 1/$pval;
													$transformPVal = log($transformPVal) if $logTransformPVal eq 'yes';
													my $score = $scoreFactor*$postionFactor*$transformPVal;
													$motifScoreHsh{$motif} += $score;
												}
											}
										}

										#---collect the score
										$strndScoreHsh{$strnd} += $motifScoreHsh{$motif};
					
										#---if mustValid but not invalid
										if ($motifScoreHsh{$motif} == 0 and $motifInfoHsh_ref->{$motif}{'mustValid'} eq 'yes') {
											$strndScoreHsh{$strnd} = 0;
											last;
										}
									}
								}
			
								if ($strndScoreHsh{'+'} > 0 or $strndScoreHsh{'-'} > 0) {
									$strndScoreHsh{$_} = sprintf "%.5f", $strndScoreHsh{$_} foreach (qw/+ -/);
									$cntgPAScoreAry_ref->[$i] = join ',', ($strndScoreHsh{'+'}, $strndScoreHsh{'-'});
									#print TMPLOG join "\t", ($strndScoreHsh{'+'}, $strndScoreHsh{'-'}."\n");
								}
							}
						
							store($cntgPAScoreAry_ref, $genomeWidePolyASitePlsPathHsh_ref->{$cntg});
						
						}
					}
				,($cntgAry_ref)
			);
		}

		#---wait until all threads are finished
		&checkRunningThreadAndWaitToJoin('yes', 1);#->490

	} else {
		
		&reportStatus("TSS Score storable found. Skip calculating", 0, "\n");#->2322
	
	}
	
	return ($genomeWidePolyASitePlsPathHsh_ref);
}
sub printCMDLogOrFinishMessage {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|626
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|69, 11_finishingTasks|300
#	secondaryAppearInSection: >none
#	input: $CMDLogOrFinishMessage
#	output: none
#	toCall: &printCMDLogOrFinishMessage($CMDLogOrFinishMessage);
#	calledInLine: 75, 305
#....................................................................................................................................................#

	my ($CMDLogOrFinishMessage) = @_;
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $absoluteScriptPath = abs_path($0);
		my $dirPath = dirname(rel2abs($absoluteScriptPath));
		my ($scriptName, $callScriptPath, $scriptSuffix) = fileparse($absoluteScriptPath, qr/\.[^.]*/);
		open (CMDLOG, ">>$dirPath/$scriptName.cmd.log.txt"); #---append the CMD log file
		print CMDLOG "[".&currentTime()."]\t"."$dirPath/$scriptName$scriptSuffix ".(join " ", @ARGV)."\n";#->626
		close CMDLOG;
		print "\n=========================================================================\n";
		print "[".&currentTime()."] starts running ...... \n";#->626
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		print "\n=========================================================================\n";
		print "[".&currentTime()."] finished running .......\n";#->626
		print "=========================================================================\n\n";
	}
}
sub printPAScoreWiggle {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: checkRunningThreadAndWaitToJoin|490, printWiggleSingleTrackFromCntgCovPlsPathHsh|2093, reportStatus|2322
#	appearInSub: >none
#	primaryAppearInSection: 9_predictTSS|271
#	secondaryAppearInSection: >none
#	input: $PAScoreWigDir, $genomeWidePolyASitePlsPathHsh_ref
#	output: 
#	toCall: &printPAScoreWiggle($PAScoreWigDir, $genomeWidePolyASitePlsPathHsh_ref);
#	calledInLine: 277
#....................................................................................................................................................#
	my ($PAScoreWigDir, $genomeWidePolyASitePlsPathHsh_ref) = @_;
	
	my %tmpStrndInfoHsh = ();
	$tmpStrndInfoHsh{'+'}{'wigPath'} = "$PAScoreWigDir/genomeWidePAScore.plus.wig";
	$tmpStrndInfoHsh{'+'}{'aryIndex'} = 0;
	$tmpStrndInfoHsh{'-'}{'wigPath'} = "$PAScoreWigDir/genomeWidePAScore.minus.wig";
	$tmpStrndInfoHsh{'-'}{'aryIndex'} = 1;
	
	&reportStatus("Printing TSS Score Wiggle", 0, "\n");#->2322

	foreach my $strnd (keys %tmpStrndInfoHsh) {
		my $wigPath = $tmpStrndInfoHsh{$strnd}{'wigPath'};
		my $aryIndex = $tmpStrndInfoHsh{$strnd}{'aryIndex'};
		my $gzip = 'no';
		my $cntgCovPlsPathHsh_ref = $genomeWidePolyASitePlsPathHsh_ref;
		if (not -s $wigPath) {
			threads->create(\&printWiggleSingleTrackFromCntgCovPlsPathHsh, ($cntgCovPlsPathHsh_ref, $aryIndex, $wigPath, $gzip));#->2093
		}
	}

	&checkRunningThreadAndWaitToJoin('yes', 1);#->490

	return ();
}
sub printWiggleSingleTrackFromCntgCovPlsPathHsh {
#....................................................................................................................................................#
#	subroutineCategory: wiggle
#	dependOnSub: >none
#	appearInSub: printPAScoreWiggle|2057
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 9_predictTSS|271
#	input: $aryIndex, $cntgCovPlsPathHsh_ref, $gzip, $wigPath
#	output: none
#	toCall: &printWiggleSingleTrackFromCntgCovPlsPathHsh($cntgCovPlsPathHsh_ref, $aryIndex, $wigPath, $gzip);
#	calledInLine: 2084
#....................................................................................................................................................#
	
	my ($cntgCovPlsPathHsh_ref, $aryIndex, $wigPath, $gzip) = @_;
	
	open (WIGGLE, ">", $wigPath);
	
	foreach my $cntg (sort keys %{$cntgCovPlsPathHsh_ref}) {

		print WIGGLE "variableStep chrom=$cntg span=1\n";

		my $cntgCovPlsPath = "$cntgCovPlsPathHsh_ref->{$cntg}";
 		system ("gzip -df $cntgCovPlsPath.gz") if (-s "$cntgCovPlsPath.gz" and $gzip eq 'yes');
		my $cntgCovAry_ref = retrieve($cntgCovPlsPath);
		system ("gzip -f $cntgCovPlsPath") if (-s $cntgCovPlsPath and $gzip eq 'yes');
		for my $i (0..$#{$cntgCovAry_ref}) {
			if ($cntgCovAry_ref->[$i]) {
				my @tmpCovAry = split /,/, $cntgCovAry_ref->[$i];
				my $cov = $tmpCovAry[$aryIndex];
				if ($cov > 0) {
					my $pos = $i + 1;
					print WIGGLE join '', ((join "\t", ($pos, $cov)), "\n");
				}
			}
		}
	}
	close WIGGLE;
}
sub readGFF_oneRNAPerGene {
#....................................................................................................................................................#
#	subroutineCategory: general, gff
#	dependOnSub: currentTime|626
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|191
#	secondaryAppearInSection: >none
#	input: $gffPath
#	output: $geneInfoHsh_ref
#	toCall: my ($geneInfoHsh_ref) = &readGFF_oneRNAPerGene($gffPath);
#	calledInLine: 203
#....................................................................................................................................................#

	my ($gffPath) = @_;

	my $geneInfoHsh_ref = {};
	
	#---read the gff
	my $geneByRNAHsh_ref = {};

	open (GFF, $gffPath);
	print "[".&currentTime()."] Reading: $gffPath\n";#->626
	while (my $theLine = <GFF>) {

		chomp $theLine;
		
		last if $theLine =~ m/^##FASTA/;
		
		if ($theLine !~ m/^\#|^\@/ and $theLine !~ m/\tsupercontig\t/) {

			my ($seq, undef, $geneCategory, $featureStart, $featureEnd, undef, $geneStrd, undef, $dscrptns) = split (/\t/, $theLine);
			
			#----assigne all non -/+ will be treated as plus
			$geneStrd = "+" if (($geneStrd ne "-") and ($geneStrd ne "+"));
			
			my @dscrptnsSplt = split /;/, $dscrptns;
			my ($unqID, $parent);
			my $geneName = "unknown";
			foreach my $theDscptn (@dscrptnsSplt) {
				if ($theDscptn =~ m/^ID=/) {$unqID = substr ($theDscptn, index ($theDscptn, "=")+1);}
				if ($theDscptn =~ m/^Parent=/) {$parent = substr ($theDscptn, index ($theDscptn, "=")+1);}
				if ($theDscptn =~ m/^description=/) {$geneName = substr ($theDscptn, index ($theDscptn, "=")+1);}
			}

			if ($geneCategory eq "gene") {#---gene
				
				my $geneID = $unqID;
				
				$geneInfoHsh_ref->{$geneID}{'strnd'} = $geneStrd;
				$geneInfoHsh_ref->{$geneID}{'cntg'} = $seq;
				$geneInfoHsh_ref->{$geneID}{'description'} = uri_unescape($geneName);
				$geneInfoHsh_ref->{$geneID}{'description'} =~ s/\+/ /g;
				@{$geneInfoHsh_ref->{$geneID}{'geneRng'}} = ($featureStart, $featureEnd);

			} elsif ($geneCategory eq "CDS") {#---Only for coding genes
				
				my $RNAID = $parent;
				my $geneID = $geneByRNAHsh_ref->{$RNAID};
				push @{$geneInfoHsh_ref->{$geneID}{'CDSRng'}}, ($featureStart, $featureEnd);
				
			} elsif ($geneCategory eq "exon") {#---exon, may be exons of alternative transcripts, wiull sort out later
				my $RNAID = $parent;
				my $geneID = $geneByRNAHsh_ref->{$RNAID};
				push @{$geneInfoHsh_ref->{$geneID}{'exonRng'}}, ($featureStart, $featureEnd);
				
			} else {#---can be tRNA, rRNA, mRNA, repRNA, ncRNA
				my $RNAID = $unqID;
				my $geneID = $parent;
				$geneByRNAHsh_ref->{$RNAID} = $geneID;
				$geneInfoHsh_ref->{$geneID}{'ctgry'} = $geneCategory;
				@{$geneInfoHsh_ref->{$geneID}{'RNARng'}} = ($featureStart, $featureEnd);
				$geneInfoHsh_ref->{$geneID}{'RNAID'} = $RNAID;
			}
		}#---end of if (($theLine !~ m/^\#|^\@/) and ($theLine !~ m/\tsupercontig\t/)) {
	}#---end of while (my $theLine = <INFILE>)
	close GFF;
	
	#---get the UTR if any
	my $minUTRLength = 10;
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {

		if (exists $geneInfoHsh_ref->{$geneID}{'CDSRng'}) {
			my $exonMin = min(@{$geneInfoHsh_ref->{$geneID}{'exonRng'}});
			my $exonMax = max(@{$geneInfoHsh_ref->{$geneID}{'exonRng'}});
			my $CDSMin = min(@{$geneInfoHsh_ref->{$geneID}{'CDSRng'}});
			my $CDSMax = max(@{$geneInfoHsh_ref->{$geneID}{'CDSRng'}});

			if ($geneInfoHsh_ref->{$geneID}{'strnd'} eq '+') {
				@{$geneInfoHsh_ref->{$geneID}{'UTR5Rng'}} = ($exonMin, $CDSMin-1) if ($CDSMin-$exonMin > $minUTRLength);
				@{$geneInfoHsh_ref->{$geneID}{'UTR3Rng'}} = ($CDSMax+1, $exonMax) if ($exonMax-$CDSMax > $minUTRLength);
			} else {
				@{$geneInfoHsh_ref->{$geneID}{'UTR3Rng'}} = ($exonMin, $CDSMin-1) if ($CDSMin-$exonMin > $minUTRLength);
				@{$geneInfoHsh_ref->{$geneID}{'UTR5Rng'}} = ($CDSMax+1, $exonMax) if ($exonMax-$CDSMax > $minUTRLength);
			}
		}
	}
	
	return ($geneInfoHsh_ref);
}
sub readMultiFasta {
#....................................................................................................................................................#
#	subroutineCategory: fasta, general
#	dependOnSub: reportStatus|2322
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|191
#	secondaryAppearInSection: >none
#	input: $fastaPath
#	output: $fastaHsh_ref
#	toCall: my ($fastaHsh_ref) = &readMultiFasta($fastaPath);
#	calledInLine: 198
#....................................................................................................................................................#

	my ($fastaPath) = @_;

	my ($seq, $seqName);
	my $fastaHsh_ref = {};
	my $i = 0;

	&reportStatus("Reading: $fastaPath", 0, "\n");#->2322
	
	open (INFILE, $fastaPath);
	chomp (my $curntLine = <INFILE>); #get the first line
	while (my $nextLine = <INFILE>) {
		chomp $nextLine;
		
		#---Only two types of line in current line, the header or seq
		if ($curntLine =~ m/^>/) {#-- header line
			my @theLineSplt = split (/\|/, $curntLine);
			$seqName = $theLineSplt[0]; #---get the first tag
			$seqName =~ s/ //g; #---remove space
			$seqName =~ s/>//g; #---remove space
		} else {#--seq line
			$seq = $seq.$curntLine;
		}
		
		#---check if next line has a > or that's the end of file
		if ($nextLine =~ m/^>/) {
			$seq =~ tr/a-z/A-Z/;
			$fastaHsh_ref->{$seqName} = $seq;
			$seq = "";

			#---ad hoc limit
			#my $cntgNum = keys %{$fastaHsh_ref}; last if $cntgNum > 100;

		} elsif (eof(INFILE)) {#---this is the last line
			$seq =~ tr/a-z/A-Z/;
			$seq = $seq.$nextLine;
			$fastaHsh_ref->{$seqName} = $seq;
		}
		
		#---next line becomes current line
		$curntLine = $nextLine;
	}

	close INFILE;
	return ($fastaHsh_ref);
}
sub readParameters {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|69
#	secondaryAppearInSection: >none
#	input: none
#	output: $fastaPath, $gffPath, $mRNABasedPolyAInfoHshPlsPath, $outDir
#	toCall: my ($mRNABasedPolyAInfoHshPlsPath, $fastaPath, $gffPath, $outDir) = &readParameters();
#	calledInLine: 78
#....................................................................................................................................................#
	
	my ($mRNABasedPolyAInfoHshPlsPath, $fastaPath, $gffPath, $outDir);

	my $dirPath = dirname(rel2abs($0));
	$outDir = "$dirPath/TSSMotifAnalyzer/";

	GetOptions 	("mRNABasedPolyAInfoHshPlsPath=s" => \$mRNABasedPolyAInfoHshPlsPath,
				 "fastaPath=s"  => \$fastaPath,
				 "gffPath=s"  => \$gffPath,
				 "outDir:s"  => \$outDir)

	or die	("Error in command line arguments\n");
	
	#---check file
	foreach my $fileToCheck ($mRNABasedPolyAInfoHshPlsPath, $fastaPath, $gffPath) {
		die "Can't read $fileToCheck" if (not -s $fileToCheck and not -s "$fileToCheck.gz");
	}

	system "mkdir -p -m 777 $outDir/";

	return($mRNABasedPolyAInfoHshPlsPath, $fastaPath, $gffPath, $outDir);
}
sub reportStatus {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|626
#	appearInSub: calculateBackgroundNucleotideFrequency|374, calculateBaseCompositionInAlignments|420, checkGeneInfo|464, checkRunningThreadAndWaitToJoin|490, createEmptyStorableForGenowideTSSPredictionData|518, createMotifHitHshStorable|561, generateShuffleSeq|763, getCoverageOfItemRngType_multiStrand|852, getCtgryGeneInfo|957, getPAScoreAroundmRNAReferencePoint|1113, getPAScoreInExonAndTSS|1207, getSequenceAroundmRNAReferencePoint|1267, getmRNAReferencePoints|1374, plotBaseCompositionAroundmRNAReferencePoint|1572, plotCovAlongCDSBothStrnd|1628, predictGenomeWidePolyASite|1903, printPAScoreWiggle|2057, readMultiFasta|2230, scanMotifAroundSiteWithMAST|2395, scanMotifWholeGenomeWithMAST|2454
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 10_analyzePAScore|282, 4_processInputData|191, 5_retrieveSequenceSurroundingPredefinedTSS|213, 6_scanMotifOccurence|235, 9_predictTSS|271
#	input: $lineEnd, $message, $numTrailingSpace
#	output: 
#	toCall: &reportStatus($message, $numTrailingSpace, $lineEnd);
#	calledInLine: 392, 397, 452, 478, 486, 513, 545, 578, 784, 793, 880, 922, 951, 976, 989, 1148, 1225, 1295, 1325, 1387, 1413, 1597, 1660, 1920, 1939, 1950, 2017, 2076, 2248, 2420, 2473, 2476, 2485, 2488, 2491
#....................................................................................................................................................#
	my ($message, $numTrailingSpace, $lineEnd) = @_;

	my $trailingSpaces = '';
	$trailingSpaces .= " " for (1..$numTrailingSpace);
	
	print "[".&currentTime()."] ".$message.$trailingSpaces.$lineEnd;#->626

	return ();
}
sub reverseComplementRefFasta {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|191
#	secondaryAppearInSection: >none
#	input: $fastaHsh_ref, $resultFastaDir
#	output: $fastaWithRevComHsh_ref, $revComFastaPath
#	toCall: my ($fastaWithRevComHsh_ref, $revComFastaPath) = &reverseComplementRefFasta($fastaHsh_ref, $resultFastaDir);
#	calledInLine: 199
#....................................................................................................................................................#
	my ($fastaHsh_ref, $resultFastaDir) = @_;
	
	my $fastaWithRevComHsh_ref = {};
	my $revComFastaPath = "$resultFastaDir/ref.with.rev.com.fasta";
	open FASTA, ">", $revComFastaPath;
	foreach my $seqName (sort keys %{$fastaHsh_ref}) {
		my $revComSeq = reverse $fastaHsh_ref->{$seqName};
		$revComSeq =~ tr/ACGTacgt/TGCAtgca/;
		print FASTA ">$seqName\+\n";
		print FASTA "$fastaHsh_ref->{$seqName}\n";
		print FASTA ">$seqName\-\n";
		print FASTA "$revComSeq\n";
		$fastaWithRevComHsh_ref->{$seqName}{'+'} = $fastaHsh_ref->{$seqName};
		$fastaWithRevComHsh_ref->{$seqName}{'-'} = $revComSeq;
	}
	close FASTA;

	return ($fastaWithRevComHsh_ref, $revComFastaPath);
}
sub runMAST {
#....................................................................................................................................................#
#	subroutineCategory: thridPartyApp
#	dependOnSub: >none
#	appearInSub: scanMotifAroundSiteWithMAST|2395, scanMotifWholeGenomeWithMAST|2454
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_scanMotifOccurence|235
#	input: $bfilePath, $extraOption, $fastaPath, $mastOutDir, $motifFilePath
#	output: $mastHitLog
#	toCall: my ($mastHitLog) = &runMAST($fastaPath, $motifFilePath, $bfilePath, $mastOutDir, $extraOption);
#	calledInLine: 2430, 2441, 2486
#....................................................................................................................................................#
	my ($fastaPath, $motifFilePath, $bfilePath, $mastOutDir, $extraOption) = @_;
	
	my $mastHitLog = "$mastOutDir/mast.hit.txt";
	my $mastErrorLog = "$mastOutDir/mast.error.txt";
	my $mastCmd = "mast $motifFilePath $fastaPath -oc $mastOutDir -ev 1e+10 -mt 10 -hit_list $extraOption >$mastHitLog 2>$mastErrorLog";
	system ("$mastCmd");

	return ($mastHitLog);
}
sub scanMotifAroundSiteWithMAST {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: getSingleMotifMASTLogPostionalData|1337, reportStatus|2322, runMAST|2374
#	appearInSub: >none
#	primaryAppearInSection: 6_scanMotifOccurence|235
#	secondaryAppearInSection: >none
#	input: $bkgdNtFreqHsh_ref, $mastRunDir, $motifFilePathHsh_ref, $predictMotifInfoHsh_ref, $resultStorableDir, $seqAroundSiteInfoHsh_ref, $shuffleSeqAroundSiteInfoHsh_ref
#	output: $mastAroundSiteResultHsh_ref
#	toCall: my ($mastAroundSiteResultHsh_ref) = &scanMotifAroundSiteWithMAST($seqAroundSiteInfoHsh_ref, $shuffleSeqAroundSiteInfoHsh_ref, $motifFilePathHsh_ref, $bkgdNtFreqHsh_ref, $mastRunDir, $predictMotifInfoHsh_ref, $resultStorableDir);
#	calledInLine: 242
#....................................................................................................................................................#
	my ($seqAroundSiteInfoHsh_ref, $shuffleSeqAroundSiteInfoHsh_ref, $motifFilePathHsh_ref, $bkgdNtFreqHsh_ref, $mastRunDir, $predictMotifInfoHsh_ref, $resultStorableDir) = @_;

	my $mastAroundSiteResultHshPlsPath = "$resultStorableDir/mastAroundSiteResultHsh.pls";
	my $mastAroundSiteResultHsh_ref = {};
	
		#my $extraOption = ' -best '; #---extrat option in MAST
		#my $extraOption = ' -comp '; #---extrat option in MAST
	my $extraOption = ' -norc ';#---extrat option in MAST
	#foreach my $siteType (keys %{$seqAroundSiteInfoHsh_ref}) {
	foreach my $siteType ('polyA_tail') {#----do the mRNA_TSS onl,
		foreach my $motif (keys %{$predictMotifInfoHsh_ref}) {
			my $motifFilePath = $motifFilePathHsh_ref->{$motif};
			my $maxHitPVal = $predictMotifInfoHsh_ref->{$motif}{'maxPValDefineBound'};
	
			&reportStatus("Running mast for $motif around $siteType", 0, "\n");#->2322
			
			my $maxPos = $seqAroundSiteInfoHsh_ref->{$siteType}{'length'};
			my $totalSeqNum = $seqAroundSiteInfoHsh_ref->{$siteType}{'totalSeqNum'};
	
			{#---on query seq
				my $fastaPath = $seqAroundSiteInfoHsh_ref->{$siteType}{'fastaPath'};
				my $subMastOutDir = "$mastRunDir/$motif/query/$siteType/";
				system ("mkdir -pm 777 $subMastOutDir");
				my $bfilePath = $bkgdNtFreqHsh_ref->{$siteType}{'full'}{'bfilePath'};
				my $mastHitLog = &runMAST($fastaPath, $motifFilePath, $bfilePath, $subMastOutDir, $extraOption);#->2374
				my ($motifPctHsh_ref, $tmpHitBySeqHsh_ref) = &getSingleMotifMASTLogPostionalData($mastHitLog, $maxHitPVal, $maxPos, $totalSeqNum);#->1337
				$mastAroundSiteResultHsh_ref->{$siteType}{$motif}{'query'}{'tmpHitBySeqHsh_ref'} = $tmpHitBySeqHsh_ref;#---take the has out of the lexcial scope
				$mastAroundSiteResultHsh_ref->{$siteType}{$motif}{'query'}{'motifPctHsh_ref'} = $motifPctHsh_ref;#---take the has out of the lexcial scope
			}
	
			{#---on shuffle seq
				my $fastaPath = $shuffleSeqAroundSiteInfoHsh_ref->{$siteType}{'full'}{'fastaPath'};
				my $subMastOutDir = "$mastRunDir/$motif/shuffle/$siteType/";
				system ("mkdir -pm 777 $subMastOutDir");
				my $bfilePath = $bkgdNtFreqHsh_ref->{$siteType}{'full'}{'bfilePath'};
				my $mastHitLog = &runMAST($fastaPath, $motifFilePath, $bfilePath, $subMastOutDir, $extraOption);#->2374
				my ($motifPctHsh_ref, $tmpHitBySeqHsh_ref) = &getSingleMotifMASTLogPostionalData($mastHitLog, $maxHitPVal, $maxPos, $totalSeqNum);#->1337
				$mastAroundSiteResultHsh_ref->{$siteType}{$motif}{'shuffle'}{'tmpHitBySeqHsh_ref'} = $tmpHitBySeqHsh_ref;#---take the has out of the lexcial scope
				$mastAroundSiteResultHsh_ref->{$siteType}{$motif}{'shuffle'}{'motifPctHsh_ref'} = $motifPctHsh_ref;#---take the has out of the lexcial scope
			}
		}
	}
	
	store($mastAroundSiteResultHsh_ref, $mastAroundSiteResultHshPlsPath);
	
	return ($mastAroundSiteResultHsh_ref);
}
sub scanMotifWholeGenomeWithMAST {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: createMotifHitHshStorable|561, generateMASTBackgroundFile|718, getMastGenomeBothStrandHit|994, reportStatus|2322, runMAST|2374, storeMotifHitToHshStorable|2501
#	appearInSub: >none
#	primaryAppearInSection: 6_scanMotifOccurence|235
#	secondaryAppearInSection: >none
#	input: $fastaHsh_ref, $fastaLengthHsh_ref, $mastBackgroundDir, $mastRunDir, $maxPolymerSize, $motifFilePathHsh_ref, $resultStorableDir, $revComFastaPath
#	output: $cntgMotifHitPlsPathHsh_ref
#	toCall: my ($cntgMotifHitPlsPathHsh_ref) = &scanMotifWholeGenomeWithMAST($revComFastaPath, $fastaHsh_ref, $motifFilePathHsh_ref, $mastRunDir, $maxPolymerSize, $mastBackgroundDir, $fastaLengthHsh_ref, $resultStorableDir);
#	calledInLine: 240
#....................................................................................................................................................#
	my ($revComFastaPath, $fastaHsh_ref, $motifFilePathHsh_ref, $mastRunDir, $maxPolymerSize, $mastBackgroundDir, $fastaLengthHsh_ref, $resultStorableDir) = @_;
	
	my ($cntgMotifHitIdxPlsPath, $cntgMotifHitPlsPathHsh_ref, $allCntgMotifHitPlsExist) = &createMotifHitHshStorable($fastaHsh_ref, $resultStorableDir);#->561

	#---will do the prediction only when cntgMotifHitIdxPls doesnt exist, checked in &createMotifHitHshStorable
	if ($allCntgMotifHitPlsExist eq 'no') {
		my $bfilePath = "$mastBackgroundDir/full.genome.freq.txt";
		if (not -s $bfilePath) {
			&reportStatus("Generating full genome background file", 0, "\n");#->2322
			my (undef, undef) = &generateMASTBackgroundFile($fastaHsh_ref, $bfilePath, $maxPolymerSize);#->718
		} else {
			&reportStatus("Full genome background file found", 0, "\n");#->443	#->2322
		}

		foreach my $motif (keys %{$motifFilePathHsh_ref}) {
			my $motifFilePath = $motifFilePathHsh_ref->{$motif};
			my $subMastOutDir = "$mastRunDir/$motif/query/fullGenome/";
			system ("mkdir -pm 777 $subMastOutDir");
			my $extraOption = ' -norc ';

			&reportStatus("Running mast for $motif on full genome", 0, "\n");#->2322
			my ($mastHitLog) = &runMAST($revComFastaPath, $motifFilePath, $bfilePath, $subMastOutDir, $extraOption);#->2374

			&reportStatus("Getting mast results for $motif", 0, "\n");#->2322
			my ($allHitBySeqHsh_ref) = &getMastGenomeBothStrandHit($mastHitLog, $fastaLengthHsh_ref);#->994

			&reportStatus("Storing mast results for $motif", 0, "\n");#->2322
			&storeMotifHitToHshStorable($cntgMotifHitPlsPathHsh_ref, $allHitBySeqHsh_ref, $motif);#->2501

			$allHitBySeqHsh_ref = {};
		}
	}

	return ($cntgMotifHitPlsPathHsh_ref);
}
sub storeMotifHitToHshStorable {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: scanMotifWholeGenomeWithMAST|2454
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_scanMotifOccurence|235
#	input: $allHitBySeqHsh_ref, $cntgMotifHitPlsPathHsh_ref, $motif
#	output: 
#	toCall: &storeMotifHitToHshStorable($cntgMotifHitPlsPathHsh_ref, $allHitBySeqHsh_ref, $motif);
#	calledInLine: 2492
#....................................................................................................................................................#
	my ($cntgMotifHitPlsPathHsh_ref, $allHitBySeqHsh_ref, $motif) = @_;
	
	foreach my $cntg (keys %{$allHitBySeqHsh_ref}) {
		my $cntgMotifHitHsh_ref = retrieve($cntgMotifHitPlsPathHsh_ref->{$cntg});
		foreach my $strnd (keys %{$allHitBySeqHsh_ref->{$cntg}}) {
			foreach my $pos (keys %{$allHitBySeqHsh_ref->{$cntg}{$strnd}}) {
				$cntgMotifHitHsh_ref->{$strnd}{$pos}{$motif} = $allHitBySeqHsh_ref->{$cntg}{$strnd}{$pos};
			}
		}
		store($cntgMotifHitHsh_ref, $cntgMotifHitPlsPathHsh_ref->{$cntg});
	}
	return ();
}

exit;
