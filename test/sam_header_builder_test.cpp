#include <sstream>

#include "sam/sam_header_builder.h"
#include "sam/sam_reader.h"

#include <boost/test/unit_test.hpp>
#include <stdexcept>

using namespace std;
using namespace gamgee;

BOOST_AUTO_TEST_CASE( add_header_line ) {
  auto reader = SingleSamReader{"testdata/test_simple.bam"};
  auto original_header = reader.header();
  auto original_header_text = original_header.header_text();

  // Replace header line.
  auto builder = SamHeaderBuilder{original_header};
  auto header_line = SamHeaderLine();
  header_line.version = "0.0";
  header_line.sorting_order = SamHeaderLine::SortingOrder::UNSORTED;
  auto modified_header = builder.set_header_line(header_line).build();

  // Verify that set_header_line() worked.
  std::stringstream header_stream;
  header_stream << modified_header.header_line();
  BOOST_CHECK_EQUAL( header_stream.str(), "@HD\tVN:0.0\tSO:unsorted" );
  // Verify that the original header is unaffected
  BOOST_CHECK_EQUAL( original_header.header_text(), original_header_text );
  // Verify that other fields were inherited from the original header.
  BOOST_CHECK( modified_header.read_groups() == original_header.read_groups());
  BOOST_CHECK( modified_header.programs() == original_header.programs());
  BOOST_CHECK( modified_header.comments() == original_header.comments());
  BOOST_CHECK_EQUAL( modified_header.n_sequences(), original_header.n_sequences());
  BOOST_CHECK_EQUAL( modified_header.sequence_name(0), original_header.sequence_name(0));
}

BOOST_AUTO_TEST_CASE( add_header_program ) {
  auto reader = SingleSamReader{"testdata/test_simple.bam"};
  auto original_header = reader.header();
  auto original_header_text = original_header.header_text();

  // Replace programs field.
  auto builder = SamHeaderBuilder{original_header};
  auto programs = std::vector<Program>();
  programs.push_back(Program());
  programs[0].id = "1";
  programs[0].command_line = "First command.";
  auto program_2 = Program();
  program_2.id = "2";
  program_2.command_line = "Second command.";
  auto modified_header = builder.set_programs(programs).append_program(program_2).build();

  // Verify that set_programs()/append_program() worked.
  BOOST_CHECK_EQUAL( modified_header.programs().size(), 2);
  BOOST_CHECK_EQUAL( modified_header.programs()[0].command_line, "First command." );
  BOOST_CHECK_EQUAL( modified_header.programs()[1].command_line, "Second command." );
  // Verify that the original header is unaffected
  BOOST_CHECK_EQUAL( original_header.header_text(), original_header_text );
  // Verify that other fields were inherited from the original header.
  BOOST_CHECK( modified_header.header_line() == original_header.header_line());
  BOOST_CHECK( modified_header.read_groups() == original_header.read_groups());
  BOOST_CHECK( modified_header.comments() == original_header.comments());
  BOOST_CHECK_EQUAL( modified_header.n_sequences(), original_header.n_sequences());
  BOOST_CHECK_EQUAL( modified_header.sequence_name(0), original_header.sequence_name(0));
}

BOOST_AUTO_TEST_CASE( add_header_read_group ) {
  auto reader = SingleSamReader{"testdata/test_simple.bam"};
  auto original_header = reader.header();
  auto original_header_text = original_header.header_text();

  // Replace read groups.
  auto builder = SamHeaderBuilder{original_header};
  auto read_groups = std::vector<ReadGroup>();
  read_groups.push_back(ReadGroup());
  read_groups[0].id = "1";
  read_groups[0].description = "First read group.";
  auto read_group_2 = ReadGroup();
  read_group_2.id = "2";
  read_group_2.description = "Second read group.";
  auto modified_header = builder.set_read_groups(read_groups).append_read_group(read_group_2).build();

  // Verify that set_read_groups()/append_read_group() worked.
  BOOST_CHECK_EQUAL( modified_header.read_groups().size(), 2);
  BOOST_CHECK_EQUAL( modified_header.read_groups()[0].description, "First read group." );
  BOOST_CHECK_EQUAL( modified_header.read_groups()[1].description, "Second read group." );
  // Verify that the original header is unaffected
  BOOST_CHECK_EQUAL( original_header.header_text(), original_header_text );
  // Verify that other fields were inherited from the original header.
  BOOST_CHECK( modified_header.header_line() == original_header.header_line());
  BOOST_CHECK( modified_header.programs() == original_header.programs());
  BOOST_CHECK( modified_header.comments() == original_header.comments());
  BOOST_CHECK_EQUAL( modified_header.n_sequences(), original_header.n_sequences());
  BOOST_CHECK_EQUAL( modified_header.sequence_name(0), original_header.sequence_name(0));
}

BOOST_AUTO_TEST_CASE( add_header_comment ) {
  auto reader = SingleSamReader{"testdata/test_simple.bam"};
  auto original_header = reader.header();
  auto original_header_text = original_header.header_text();

  // Replace comments.
  auto builder = SamHeaderBuilder{original_header};
  auto comments = std::vector<SamHeaderComment>();
  comments.push_back(SamHeaderComment());
  comments[0].comment = "First comment.";
  auto comment_2 = SamHeaderComment();
  comment_2.comment = "Second comment.";
  auto modified_header = builder.set_header_comments(comments).append_header_comment(comment_2).build();

  // Verify that set_comments()/append_comment() worked.
  BOOST_CHECK_EQUAL( modified_header.comments().size(), 2);
  BOOST_CHECK_EQUAL( modified_header.comments()[0].comment, "First comment." );
  BOOST_CHECK_EQUAL( modified_header.comments()[1].comment, "Second comment." );
  // Verify that the original header is unaffected
  BOOST_CHECK_EQUAL( original_header.header_text(), original_header_text );
  // Verify that other fields were inherited from the original header.
  BOOST_CHECK( modified_header.header_line() == original_header.header_line());
  BOOST_CHECK( modified_header.programs() == original_header.programs());
  BOOST_CHECK( modified_header.read_groups() == original_header.read_groups());
  BOOST_CHECK_EQUAL( modified_header.n_sequences(), original_header.n_sequences());
  BOOST_CHECK_EQUAL( modified_header.sequence_name(0), original_header.sequence_name(0));
}

BOOST_AUTO_TEST_CASE( add_header_seq_info ) {
  auto reader = SingleSamReader{"testdata/test_simple.bam"};
  auto original_header = reader.header();
  auto original_header_text = original_header.header_text();

  // Replace sequences info.
  auto builder = SamHeaderBuilder{original_header};
  auto seqs_info = std::vector<std::pair<std::string, uint32_t>>();
  seqs_info.push_back(std::make_pair("seq1", 1));
  seqs_info.push_back(std::make_pair("seq2", 2));
  auto modified_header = builder.set_seqs_info(seqs_info).append_seq_info("seq3", 3).build();
  // Verify that set_programs()/append_program() worked.
  BOOST_CHECK_EQUAL( modified_header.n_sequences(), 3);
  BOOST_CHECK_EQUAL( modified_header.sequence_name(0), "seq1" );
  BOOST_CHECK_EQUAL( modified_header.sequence_length(1), 2);
  // Verify that the original header is unaffected
  BOOST_CHECK_EQUAL( original_header.header_text(), original_header_text );
  // Verify that other fields were inherited from the original header.
  BOOST_CHECK( modified_header.header_line() == original_header.header_line());
  BOOST_CHECK( modified_header.read_groups() == original_header.read_groups());
  BOOST_CHECK( modified_header.programs() == original_header.programs());
  BOOST_CHECK( modified_header.comments() == original_header.comments());
}



BOOST_AUTO_TEST_CASE( build_multiple_headers ) {
  auto reader = SingleSamReader{"testdata/test_simple.bam"};
  auto original_header = reader.header();
  auto original_header_text = original_header.header_text();

  // Replace comments.
  auto builder = SamHeaderBuilder{original_header};
  auto comments = std::vector<SamHeaderComment>();
  comments.push_back(SamHeaderComment());
  comments[0].comment = "Comment.";
  auto modified_header = builder.set_header_comments(comments).build();

  BOOST_CHECK( modified_header.comments() == comments);
  // Verify that the original header is unaffected
  BOOST_CHECK_EQUAL( original_header.header_text(), original_header_text );

  // Add read group.
  auto read_group = ReadGroup();
  read_group.id = "1";
  read_group.description = "Read group.";
  auto modified_header_2 = builder.append_read_group(read_group).build();
  auto n_original_read_groups = original_header.read_groups().size();

  BOOST_CHECK( modified_header_2.comments() == comments);
  BOOST_CHECK_EQUAL( modified_header_2.read_groups().size(),
                     n_original_read_groups + 1);
  BOOST_CHECK_EQUAL( modified_header_2.read_groups()[n_original_read_groups].description,
                     "Read group." );
  // Check previous header is unchanged
  BOOST_CHECK_EQUAL( modified_header.read_groups().size(),
                     n_original_read_groups);
  // Verify that the original header is unaffected
  BOOST_CHECK_EQUAL( original_header.header_text(), original_header_text );
}
