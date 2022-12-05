use std::error::Error;
use std::fmt::Formatter;

#[derive(Debug)]
pub struct InputError(String);

impl InputError {
    pub fn new(err: &str) -> Self {
        Self(err.to_owned())
    }
}

impl std::fmt::Display for InputError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}
impl Error for InputError {}

impl From<InputError> for String {
    fn from(err: InputError) -> Self {
        err.0
    }
}

impl From<String> for InputError {
    fn from(s: String) -> Self {
        Self(s)
    }
}

impl From<&str> for InputError {
    fn from(s: &str) -> Self {
        Self::new(s)
    }
}
