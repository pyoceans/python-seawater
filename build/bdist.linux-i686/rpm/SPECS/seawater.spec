%define name seawater
%define version 1.0.2
%define unmangled_version 1.0.2
%define release 1

Summary: Seawater Libray for Python
Name: %{name}
Version: %{version}
Release: %{release}
Source0: %{name}-%{unmangled_version}.tar.gz
License: MIT
Group: Development/Libraries
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-buildroot
Prefix: %{_prefix}
BuildArch: noarch
Vendor: Filipe Fernandes <ocefpaf@gmail.com>
Url: http://ocefpaf.tiddlyspot.com/#python-seawater

%description
    This module (version 1.0.2) is a translation of the original SEAWATER 3.2 MATLAB toolkit routines for calculating the properties of sea water. They are a self contained library and are extremely easy to use.

%prep
%setup -n %{name}-%{unmangled_version}

%build
python setup.py build

%install
python setup.py install --root=$RPM_BUILD_ROOT --record=INSTALLED_FILES

%clean
rm -rf $RPM_BUILD_ROOT

%files -f INSTALLED_FILES
%defattr(-,root,root)
