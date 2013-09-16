using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Bio;
using Bio.IO.FastQ;
using System.IO;
using System.Text.RegularExpressions;

namespace Bio.IO
{
    /// <summary>
    /// An IlluminaParser parses a FastQ data source, much like FastQParser,
    /// and also parses the Illumina sequence identifiers.
    /// </summary>
    /// <remarks>
    /// The following information is based on the FastQ Wikipedia article: http://en.wikipedia.org/wiki/FASTQ_format
    /// and the CASAVA 1.8.2 User Guide pp.39-41: http://support.illumina.com/sequencing/sequencing_software/casava/documentation.ilmn
    /// 
    /// If the sequence identifier comes from an older version of Casava (before 1.8), the following
    /// fields are stored in the returned QualitativeSequence objects' metadata:
    /// <list type="table">
    ///     <listheader>
    ///         <term>Field</term>
    ///         <description>Description</description>
    ///     </listheader>
    ///     <item>
    ///         <term>Instrument</term>
    ///         <description><c>string</c> The unique instrument name</description>
    ///     </item>
    ///     <item>
    ///         <term>Lane</term>
    ///         <description><c>int</c> The flowcell lane</description>
    ///     </item>
    ///     <item>
    ///         <term>Tile</term>
    ///         <description><c>int</c> The tile number within the flowcell lane</description>
    ///     </item>
    ///     <item>
    ///         <term>X</term>
    ///         <description><c>int</c> The x-coordinate of the cluster within the tile</description>
    ///     </item>
    ///     <item>
    ///         <term>Y</term>
    ///         <description><c>int</c> The y-coordinate of the cluster within the tile</description>
    ///     </item>
    ///     <item>
    ///         <term>Index</term>
    ///         <description><c>int</c> The index number for a multi-plexed sample (0 for no indexing)</description>
    ///     </item>
    ///     <item>
    ///         <term>PairMember</term>
    ///         <description><c>int</c> The member of a pair, 1 or 2 (paired-end or mate-pair reads only)</description>
    ///     </item>
    /// </list>
    /// 
    /// If the sequence identifier comes from a later version of Casava (after 1.8), the following
    /// fields are stored in the returned QualitativeSequence objects' metadata:
    /// <list type="table">
    ///     <listheader>
    ///         <term>Field</term>
    ///         <description>Description</description>
    ///     </listheader>
    ///     <item>
    ///         <term>Instrument</term>
    ///         <description><c>string</c> The unique instrument name</description>
    ///     </item>
    ///     <item>
    ///         <term>Run</term>
    ///         <description><c>int</c> The run ID</description>
    ///     </item>
    ///     <item>
    ///         <term>FlowCell</term>
    ///         <description><c>string</c> The flowcell ID</description>
    ///     </item>
    ///     <item>
    ///         <term>Lane</term>
    ///         <description><c>int</c> The flowcell name</description>
    ///     </item>
    ///     <item>
    ///         <term>Tile</term>
    ///         <description><c>int</c> The tile number within the flowcell lane</description>
    ///     </item>
    ///     <item>
    ///         <term>X</term>
    ///         <description><c>int</c> The x-coordinate of the cluster within the tile</description>
    ///     </item>
    ///     <item>
    ///         <term>Y</term>
    ///         <description><c>int</c> The y-coordinate of the cluster within the tile</description>
    ///     </item>
    ///     <item>
    ///         <term>PairMember</term>
    ///         <description><c>int</c> The member of a pair, 1 or 2 (paired-end or mate-pair reads only)</description>
    ///     </item>
    ///     <item>
    ///         <term>IsFiltered</term>
    ///         <description><c>bool</c> Y if the read fails filter (read is bad), N otherwise</description>
    ///     </item>
    ///     <item>
    ///         <term>ControlBits</term>
    ///         <description><c>int</c> 0 when none of the control bits are on, otherwise it is an even number</description>
    ///     </item>
    ///     <item>
    ///         <term>IndexSequence</term>
    ///         <description><c>string</c> Y if the read fails filter (read is bad), N otherwise</description>
    ///     </item>
    /// </list>
    /// </remarks>
    public class IlluminaParser : FastQParser
    {
        /**
         * Before Casava version 1.8, the ID string had the following format (sans whitespace):
         * @ Instrument : Lane : Tile : X : Y # Index / PairMember
         */
        private const string pre18Pattern = "^@([^:]+):([\\d]+):([\\d]+):([\\d]+):([\\d]+)#([\\d]+/([12])$";
        /**
         * After Casava version 1.8, the ID string has the following format:
         * @ Instrument : Run : FlowCell : Lane : Tile : X : Y    PairMember : IsFiltered : ControlBits : IndexSequence
         */
        private const string post18Pattern = "^@([a-zA-Z0-9_]+):([\\d]+):([a-zA-Z0-9]+):([\\d]+):([\\d]+):([\\d]+):([\\d]+) ([12]):([YN]):([\\d]+):([ACGT]+)$";

        private static Regex pre18Regex;
        private static Regex post18Regex;

        protected static Regex Pre18Regex
        {
            get
            {
                if (pre18Regex == null)
                {
                    pre18Regex = new Regex(pre18Pattern);
                }
                return pre18Regex;
            }
        }

        protected static Regex Post18Regex
        {
            get
            {
                if (post18Regex == null)
                {
                    post18Regex = new Regex(post18Pattern);
                }
                return post18Regex;
            }
        }

        public IlluminaParser() : base() { }
        public IlluminaParser(string filename) : base(filename) { }

        public override IEnumerable<QualitativeSequence> Parse()
        {
            using (StreamReader streamReader = new StreamReader(this.Filename))
            {
                FastQFormatType formatType = this.FormatType;
                do
                {
                    var seq = ParseOne(streamReader, formatType);
                    if (seq != null)
                        yield return ParseHeader(seq);
                }
                while (!streamReader.EndOfStream);
            }
        }

        /* Cannot override, because base is not virtual
        public override IEnumerable<QualitativeSequence> Parse(StreamReader reader)
        {
            FastQFormatType formatType = this.FormatType;
            do
            {
                var seq = ParseOne(reader, formatType);
                if (seq != null)
                    yield return ParseHeader(seq);
            }
            while (!reader.EndOfStream);
        }*/

        private QualitativeSequence ParseHeader(QualitativeSequence seq)
        {
            Match m = pre18Regex.Match(seq.ID);
            if (m.Success)
            {
                seq.Metadata["Instrument"]  = m.Captures[0].Value;
                seq.Metadata["Lane"]        = Int32.Parse(m.Captures[1].Value);
                seq.Metadata["Tile"]        = Int32.Parse(m.Captures[2].Value);
                seq.Metadata["X"]           = Int32.Parse(m.Captures[3].Value);
                seq.Metadata["Y"]           = Int32.Parse(m.Captures[4].Value);
                seq.Metadata["Index"]       = Int32.Parse(m.Captures[5].Value);
                seq.Metadata["PairMember"]  = Int32.Parse(m.Captures[6].Value);
                return seq;
            }
            else
            {
                m = post18Regex.Match(seq.ID);
                if (m.Success)
                {
                    seq.Metadata["Instrument"]      = m.Captures[0].Value;
                    seq.Metadata["Run"]             = Int32.Parse(m.Captures[1].Value);
                    seq.Metadata["FlowCell"]        = m.Captures[2].Value;
                    seq.Metadata["Lane"]            = Int32.Parse(m.Captures[3].Value);
                    seq.Metadata["Tile"]            = Int32.Parse(m.Captures[4].Value);
                    seq.Metadata["X"]               = Int32.Parse(m.Captures[5].Value);
                    seq.Metadata["Y"]               = Int32.Parse(m.Captures[6].Value);
                    seq.Metadata["PairMember"]      = Int32.Parse(m.Captures[7].Value);
                    seq.Metadata["IsFiltered"]      = m.Captures[8].Value == "Y";
                    seq.Metadata["ControlBits"]     = Int32.Parse(m.Captures[9].Value);
                    seq.Metadata["IndexSequence"]   = m.Captures[10].Value;
                    return seq;
                }
            }
            throw new FileFormatException("Sequence identifier not in Illumina format");
        }
    }


}
