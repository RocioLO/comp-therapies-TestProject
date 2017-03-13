project(org_mitk_points_registration)

mitk_create_plugin(
  EXPORT_DIRECTIVE REGISTRATION_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkQtWidgetsExt
)
