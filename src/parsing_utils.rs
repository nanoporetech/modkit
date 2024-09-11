use nom::bytes::complete::tag;
use nom::character::complete::{
    alphanumeric1, anychar, multispace1, none_of, u64 as nomu64,
};
use nom::multi::{fold_many1, separated_list0};
use nom::number::complete::float;
use nom::IResult;

pub(crate) fn consume_digit(l: &str) -> IResult<&str, u64> {
    multispace1(l).and_then(|(r, _)| nomu64(r))
}

pub(crate) fn consume_float(l: &str) -> IResult<&str, f32> {
    multispace1(l).and_then(|(r, _)| float(r))
}

pub(crate) fn consume_dot(l: &str) -> IResult<&str, &str> {
    multispace1(l).and_then(|(r, _)| tag(".")(r))
}

#[allow(dead_code)] // keeping this in case I want it later.. so I don't have to reinvent it
pub(crate) fn consume_char_from_list<'a>(
    l: &'a str,
    sep: &str,
) -> IResult<&'a str, char> {
    separated_list0(tag(sep), alphanumeric1)(l)
        .and_then(|(r, parts)| anychar(parts[0]).map(|(_, mc)| (r, mc)))
}

pub(crate) fn consume_string_from_list<'a>(
    l: &'a str,
    sep: &str,
) -> IResult<&'a str, &'a str> {
    separated_list0(tag(sep), alphanumeric1)(l).map(|(r, parts)| (r, parts[0]))
}

pub(crate) fn consume_string_spaces(l: &str) -> IResult<&str, String> {
    fold_many1(none_of("\t\r\n"), String::new, |mut acc, item| {
        acc.push(item);
        acc
    })(l)
}

pub(crate) fn consume_string(l: &str) -> IResult<&str, String> {
    fold_many1(none_of(" \t\r\n"), String::new, |mut acc, item| {
        acc.push(item);
        acc
    })(l)
}

pub(crate) fn consume_char(l: &str) -> IResult<&str, char> {
    multispace1(l).and_then(|(r, _)| anychar(r))
}
