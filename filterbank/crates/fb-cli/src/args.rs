use std::collections::HashMap;

pub struct Args {
    pub positional: Vec<String>,
    pub options: HashMap<String, String>,
    pub flags: Vec<String>,
}

impl Args {
    pub fn parse(raw: &[String], with_value: &[&str]) -> Result<Self, String> {
        let mut positional = Vec::new();
        let mut options = HashMap::new();
        let mut flags = Vec::new();
        let mut i = 0;
        while i < raw.len() {
            let a = &raw[i];
            if let Some(name) = a.strip_prefix("--").or_else(|| a.strip_prefix('-').filter(|s| s.len() == 1)) {
                let (k, inline) = match name.split_once('=') {
                    Some((k, v)) => (k.to_string(), Some(v.to_string())),
                    None => (name.to_string(), None),
                };
                if with_value.contains(&k.as_str()) {
                    let v = match inline {
                        Some(v) => v,
                        None => {
                            i += 1;
                            raw.get(i).cloned().ok_or_else(|| format!("--{k} needs a value"))?
                        }
                    };
                    options.insert(k, v);
                } else {
                    flags.push(k);
                }
            } else {
                positional.push(a.clone());
            }
            i += 1;
        }
        Ok(Self { positional, options, flags })
    }

    pub fn get(&self, k: &str) -> Option<&str> {
        self.options.get(k).map(|s| s.as_str())
    }

    pub fn flag(&self, k: &str) -> bool {
        self.flags.iter().any(|f| f == k)
    }

    pub fn num<T: std::str::FromStr>(&self, k: &str, default: T) -> Result<T, String> {
        match self.get(k) {
            Some(v) => v.parse().map_err(|_| format!("--{k}: cannot parse '{v}'")),
            None => Ok(default),
        }
    }
}
